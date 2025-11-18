#' Fit a Generalized Additive Partial-Mastery Cognitive Diagnosis Model (GaPM-CDM).
#'
#' @param data Matrix \code{(n x p)} with binary entries. Missing data should be coded as \code{NA}.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param start.par (optional) For simulation use. List of size 4 with starting model parameters (A, C, D, R).
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{A}: A matrix \code{(p x q)} of estimated contribution ('weights') parameters.
#'  \item \code{C}: An array \code{(p x tp x q)} of estimated I/B-spline coefficients \code{(tp = degree + |knots|)}.
#'  \item \code{R}: A matrix \code{(q x q)} of estimated correlations for the latent variables
#'  \item \code{llk}: Marginal log-likelihood (computed via Monte-Carlo integration) evaluated at \code{A}, \code{C}, and \code{R}.
#'  \item \code{AIC}: Akaike information criterion for the estimated model.
#'  \item \code{BIC}: Bayesian (Schwarz) information criterion for the estimated model.
#'  \item \code{cdllk.trace}: (if return.trace = T) A matrix \code{(iter x 2)} with the trace for the observed variables and latent variables log-likelihood.
#'  \item \code{ar.trace}: (if return.trace = T) A vector \code{(iter)} with the trace for the acceptance rate (for MALA and RWMH samplers).
#'  \item \code{theta.trace}: (if return.trace = T) A matrix \code{(iter x dim(theta))} with the trace for the parameter estimates.
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
gapmCDM <- function(data,q,control = list(), start.par = NULL, ...){

  control = pr_control_gaCDM(control,...)
  p = ncol(data)
  if(control$verbose) cat(" Model: Generalized Additive PM-CDM \n")
  control$sampler = match.arg(control$sampler, c("ULA","MALA","RWMH"))
  control$basis = match.arg(control$basis, c("is","bs","pwl"))
  control$algorithm = match.arg(control$algorithm, c("GD","ADAM","mixed"))
  if(is.null(control$degree) && control$basis == "pwl") control$degree = 1
  if(is.null(control$degree) && control$basis != "pwl") control$degree = 2
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)
  tp = length(control$knots) + control$degree

  if(!is.null(start.par) && !is.list(start.par) && (start.par == "random")){
    if(control$verbose) cat(" Generating random starting values for model parameters ...")
    control$prob.sparse = 0.0
    control$iden.R = T
    if(!is.null(control$seed)) set.seed(control$seed)
    pp = pr_param_gaCDM(p,q,tp,T,control)
    if(control$verbose) cat("\r Generating random starting values for model parameters ... (Done!) \n")
  }
  if(is.null(start.par) || is.list(start.par)){
    pp = pr_param_gaCDM(p,q,tp,F,control)
    if(is.list(start.par)){
      if(!all(names(start.par) %in% names(pp))) stop("Argument `start.par' must be a list with elements `A', `C', `D', `mu', and/or `R'.")
      for(ii in names(start.par)){
        if(ii == "A") if(nrow(start.par[[ii]]) != ncol(data)) stop("Matrix `start.par$A' mismatch rows with p.")
        if(ii == "C") if(ncol(start.par[[ii]]) != tp) stop("Matrix `start.par$C' mismatch cols with length(control$knots) + control$degree")
        if(ii == "D") if(ncol(start.par[[ii]]) != tp) stop("Matrix `start.par$D' mismatch cols with length(control$knots) + control$degree")
        if(ii == "mu") if(length(start.par[[ii]]) != q) stop("Lenght of `start.par$mu' mismatch with K.")
        if(ii == "R"){
          if((nrow(start.par[[ii]]) != ncol(start.par[[ii]])) || (max(nrow(start.par[[ii]]), ncol(start.par[[ii]])) != q)){
            stop("Rows or cols of `start.par$R' mismatch with K.")
          }
        }
        pp[[ii]] = start.par[[ii]]
      }
    }
  }

  if(!is.null(control$start.zn) && is.matrix(control$start.zn)){
    zn = control$start.zn
  } else if(control$start.zn == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn = mvtnorm::rmvnorm(nrow(data),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressWarnings(psych::fa(r = data, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                     missing = T))
    zn = tmp$scores
    if(control$cor.R) pp$R = cor(zn) else pp$R = cov(zn)
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  if(!is.null(control$seed)) set.seed(control$seed)
  out <- gapmCDM_fit_rcpp(Y = data[],A = pp$A[],C = pp$C[],D = pp$D[],mu = pp$mu[],R = pp$R[], Z = zn[], control = control)
  colnames(out$PI) <- rownames(out$A) <- rownames(out$C) <- colnames(data)
  colnames(out$Z) <- colnames(out$A) <- colnames(out$R) <- rownames(out$R) <- names(out$mu) <- paste0("Z",1:q)
  colnames(out$U) <- paste0("U",1:q)
  if(control$return.trace){
    colnames(out$cdllk.trace) <- c("fz","fyz")
    Anames <- paste0("A",apply(expand.grid(paste0("j",1:p),paste0("k",1:q)),1,paste,collapse = "."))
    Cnames <- paste0("C",apply(expand.grid(apply(expand.grid(paste0("j",1:p),paste0("r",1:ncol(out$C))),1,paste,collapse = "."), paste0("k",1:q)),1,paste,collapse = "."))
    Mnames <- paste0("mu",1:q)
    Rnames <- paste0("R",apply(which(lower.tri(diag(q)) == T,arr.ind = T),1,paste0,collapse = ""))
    colnames(out$theta.trace) <- c(Anames,Cnames,Mnames,Rnames)
  }
  class(out) = c("gapmCDM", "pmCDM")
  return(out)
}

#' Simulate data from a Generalized Additive Partial-Mastery Cognitive Diagnosis Model (GaPM-CDM).
#'
#' @param n Number of simulated entries.
#' @param p Number of observed (binary) variables.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param start.par (optional) List of size 4 with starting model parameters (A, C, D, R).
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{Y}: A matrix \code{(n x p)} of simulated observed (binary) variables.
#'  \item \code{Z}: A matrix \code{(n x q)} of simulated latent variables (on the continuous scale).
#'  \item \code{U}: A matrix \code{(n x q)} of simulated latent variables (on the \code{[0,1]} scale).
#'  \item \code{spM}: A matrix \code{(n x tp)} of I/B-spline basis functions (\code{(tp = degree + |knots|)}).
#'  \item \code{PI}: A matrix \code{(n x p)} of predicted probabilities.
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
gapmCDM_sim <- function(n, p, q, control = list(),
                        start.par = NULL, ...){
  control = pr_controlsim_gaCDM(control,q,...)
  control$basis = match.arg(control$basis, c("is","bs","pwl"))
  if(is.null(control$degree) && control$basis == "pwl") control$degree = 1
  if(is.null(control$degree) && control$basis != "pwl") control$degree = 2
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)
  tp = length(control$knots) + control$degree
  if(!is.null(start.par) && !is.list(start.par) && (length(start.par) != 5))
    stop("Argument `start.par' needs to be a list with elements `A', `C', `D', `mu', and `R'.")
  if(!is.null(start.par)){
    pp = start.par
    if(nrow(pp$A) != p) stop("Matrix `start.par$A' mismatch rows with p.")
    if(ncol(pp$C) != ncol(pp$D)) stop("Mismatch columns in `start.par$C' and `start.par$D'.")
    if(ncol(pp$C) != tp) stop("Matrix `start.par$C' mismatch rows with length(control$knots) + control$degree (+1 if intercept)")
    if(length(pp$mu) != nrow(pp$R)) stop("Lenght of `start.par$mu' and nrows of `start.par$R' differ.")
    if(sum(pp$A * control$Qmatrix) != p) stop("Check configuration of `start.par$A' and `control$Qmatrix'.")
  } else {
    pp = pr_param_gaCDM(p,q,tp,T,control)
  }
  if(!is.null(control$seed)) set.seed(control$seed)
  out <- gapmCDM_sim_rcpp(n,q,p,pp$A[,,drop = F],pp$C,pp$mu,pp$R,control)
  out$A <- pp$A
  out$C <- pp$C
  out$D <- pp$D
  out$mu <- pp$mu
  out$R <- pp$R
  out$posMu <- matrix(colMeans(out$Z),nrow = n, ncol = q, byrow = T)
  out$posR <- array(cov(out$Z),dim = c(q,q,n))
  rownames(out$A) <- rownames(out$C) <- rownames(out$D) <- colnames(out$Y) <- colnames(out$PI) <- paste0("Y",1:p)
  colnames(out$A) <- colnames(out$R) <- rownames(out$R) <- colnames(out$posR) <- rownames(out$posR) <- colnames(out$Z) <- names(out$mu) <- colnames(out$posMu) <- paste0("Z",1:q)
  colnames(out$U) <- paste0("U",1:q)
  class(out) = c("gapmCDM", "pmCDM")
  return(out)
}

#' Find number of latent variables (q) using cross-validation (GaPM-CDM)
#'
#' @param Ytrain Matrix with observed binary entries to train the model. Missing entries (\code{Ytest}) coded as \code{NA}.
#' @param Ytest Matrix with observed binary entries to test the model. All entries are \code{NA} but the missing in \code{Ytrain}.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{CV.error}: Cross-validation error.
#'  \item \code{AUC}: Area under the curve error.
#'  \item \code{llk.gapmCDM} Marginal log-likelihood for the gapmCDM model
#' }
#' @details Define CV.error.
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
gapmCDM_findqCV <- function(Ytrain, Ytest, q, control = list(), start.par = NULL, ...){

  control = pr_control_gaCDM(control,...)
  if(control$verbose) cat(" Model: Generalized Additive PM-CDM \n")
  control$sampler = match.arg(control$sampler, c("ULA","MALA","RWMH"))
  control$basis = match.arg(control$basis, c("is","bs","pwl"))
  control$algorithm = match.arg(control$algorithm, c("GD","ADAM","mixed"))
  if(is.null(control$degree) && control$basis == "pwl") control$degree = 1
  if(is.null(control$degree) && control$basis != "pwl") control$degree = 2
  p = ncol(data)
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)
  tp = length(control$knots) + control$degree

  if(is.list(start.par)){
    if(length(start.par) != 5) stop("Argument `start.par' needs to be a list with elements `A', `C', `D', `mu', and `R'.")
    ppD1 = start.par
    if(nrow(ppD1$A) != ncol(data)) stop("Matrix `start.par$A' mismatch rows with p.")
    if(ncol(ppD1) != ncol(ppD1$D)) stop("Mismatch columns in `start.par$C' and `start.par$D'.")
    if(ncol(ppD1$C) != tp) stop("Matrix `start.par$C' mismatch rows with length(control$knots) + control$degree (+1 if intercept)")
    if(length(ppD1$mu) != nrow(ppD1$R)) stop("Lenght of `start.par$mu' and nrows of `start.par$R' differ.")
  }
  if(!is.null(start.par) && !is.list(start.par) && (start.par == "random")){
    if(control$verbose) cat(" Generating random starting values for model parameters ...")
    control$prob.sparse = 0.0
    control$iden.R = T
    if(!is.null(control$seed)) set.seed(control$seed)
    ppD1 = pr_param_gaCDM(p,q,tp,T,control)
    if(control$verbose) cat("\r Generating random starting values for model parameters ... (Done!) \n")
  }
  if(is.null(start.par)){
    ppD1 = pr_param_gaCDM(p,q,tp,F,control)
  }

  if(!is.null(control$start.zn) && is.matrix(control$start.zn)){
    zn = control$start.zn
  } else if(control$start.zn == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn = mvtnorm::rmvnorm(nrow(Ytrain),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressWarnings(psych::fa(r = Ytrain, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                     missing = T))
    zn = tmp$scores
    if(control$cor.R) ppD1$R = cor(zn) else ppD1$R = cov(zn)
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  if(!is.null(control$seed)) set.seed(control$seed)
  if(control$verbose) cat(paste0("\n [Training data] \n"))
  fit1 = gapmCDM_fit_rcpp(Y = Ytrain[],A = ppD1$A[],C = ppD1$C[],D = ppD1$D[],mu = ppD1$mu[], R = ppD1$R[], Z = zn[], control = control)
  out1 = gapmCDM_cv_rcpp(Ytrain[], Ytest[],A = fit1$A[],C = fit1$C[],mu = fit1$mu[],R = fit1$R[], Z = fit1$Z[], control = control)
  pred = ROCR::prediction(out1$Yhat,out1$Yobs)
  aucCV = unlist(methods::slot(ROCR::performance(pred, "auc"), "y.values"))

  return(list(CV.error = out1$CV.error, AUC = aucCV, llk.train = fit1$llk,
              mod.train = fit1))
}


#' Find number of latent variables (q) using cross-validated MLLK
#'
#' @param Ytrain Matrix with observed binary entries to train the model. Missing entries (\code{Ytest}) coded as \code{NA}.
#' @param Ytest Matrix with observed binary entries to test the model. All entries are \code{NA} but the missing in \code{Ytrain}.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{mllk.test}: Test data marginal log-likelihood.
#'  \item \code{mllk.train}: Train data marginal log-likelihood.
#'  \item \code{train.mod} Fitted model on Train data (for traceability).
#' }
#' @details Define CV.error.
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
gapmCDM_mllkCV <- function(Ytrain, Ytest, q, control = list(), ...){
  p = ncol(Ytrain)
  control = pr_control_gaCDM(control,...)
  if(control$verbose) cat(" Model: Generalized Additive PM-CDM \n")
  if(control$verbose) cat(paste0("\n [ Training data ] \n\n"))
  control$sampler = match.arg(control$sampler, c("ULA","MALA","RWMH"))
  control$basis = match.arg(control$basis, c("is","bs","pwl"))
  if(is.null(control$degree) && control$basis == "pwl") control$degree = 1
  if(is.null(control$degree) && control$basis != "pwl") control$degree = 2
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)

  tp = length(control$knots) + control$degree
  ppD1 = pr_param_gaCDM(p,q,tp,F,control)

  if(!is.null(control$start.zn) && is.matrix(control$start.zn)){
    zn = control$start.zn
  } else if(control$start.zn == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn = mvtnorm::rmvnorm(nrow(Ytrain),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressMessages(suppressWarnings(psych::fa(r = Ytrain, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                                      missing = T)))
    zn = tmp$scores
    if(control$cor.R) ppD1$R = cor(zn) else ppD1$R = cov(zn)
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  if(!is.null(control$seed)) set.seed(control$seed)
  fit1 = gapmCDM_fit_rcpp(Y = Ytrain[],A = ppD1$A[],C = ppD1$C[],D = ppD1$D[],mu = ppD1$mu[],R = ppD1$R[], Z = zn[], control = control)
  if(control$return.trace){
    colnames(fit1$cdllk.trace) <- c("fz","fyz")
    Anames <- paste0("A",apply(expand.grid(paste0("j",1:p),paste0("k",1:q)),1,paste,collapse = "."))
    Cnames <- paste0("C",apply(expand.grid(apply(expand.grid(paste0("j",1:p),paste0("r",1:ncol(fit1$C))),1,paste,collapse = "."), paste0("k",1:q)),1,paste,collapse = "."))
    Mnames <- paste0("mu",1:q)
    Rnames <- paste0("R",apply(which(lower.tri(diag(q)) == T,arr.ind = T),1,paste0,collapse = ""))
    colnames(fit1$theta.trace) <- c(Anames,Cnames,Mnames,Rnames)
  }
  if(control$verbose) cat(paste0("\n\n [ Testing data ] \n\n"))

  if(!is.null(control$start.zn.test) && is.matrix(control$start.zn.test)){
    zn.test = control$start.zn.test
  } else if(control$start.zn.test == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn.test = mvtnorm::rmvnorm(nrow(Ytest),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn.test == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressMessages(suppressWarnings(psych::fa(r = Ytest, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                                      missing = T)))
    zn.test = tmp$scores
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  if(!control$cv.useFitPos){
    IS_posMu = matrix(c(fit1$mu), nrow = nrow(Ytest), ncol = q)
    IS_posR = array(fit1$R, dim = c(q,q,nrow(Ytest)))
  } else {
    IS_posMu = fit1$posMu
    IS_posR =  fit1$posR
  }

  testmllk = fy_gapmCDM_IS(Ytest[],A = fit1$A[],C = fit1$C[],mu = fit1$mu[],L = t(chol(fit1$R[])), Z = zn.test[],
                           pM = IS_posMu[], pR = IS_posR[], control = control)
  return(list(mllk.test = testmllk, mllk.train = fit1$llk, train.mod = fit1)) #
}

#' Fit an Partial-Mastery Additive Cognitive Diagnosis Model (PM-ACDM).
#'
#' @param data Matrix \code{(n x p)} with binary entries. Missing data should be coded as \code{NA}.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param start.par (optional) For simulation use. List of size 4 with starting model parameters (G, mu, R).
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{G}: A matrix \code{(p x (q+1))} of estimated parameters.
#'  \item \code{mu}: An vector \code{(q x 1)} of estimated latent variable means.
#'  \item \code{R}: A matrix \code{(q x q)} of estimated correlations for the latent variables
#'  \item \code{llk}: Marginal log-likelihood (computed via Monte-Carlo integration) evaluated at \code{G}, \code{mu}, and \code{R}.
#'  \item \code{AIC}: Akaike information criterion for the estimated model.
#'  \item \code{BIC}: Bayesian (Schwarz) information criterion for the estimated model.
#'  \item \code{cdllk.trace}: (if return.trace = T) A matrix \code{(iter x 2)} with the trace for the observed variables and latent variables log-likelihood.
#'  \item \code{ar.trace}: (if return.trace = T) A vector \code{(iter)} with the trace for the acceptance rate (for MALA and RWMH samplers).
#'  \item \code{theta.trace}: (if return.trace = T) A matrix \code{(iter x dim(theta))} with the trace for the parameter estimates.
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
apmCDM <- function(data, q, control = list(), start.par = NULL, ...){

  control = pr_control_aCDM(control,...)
  p = ncol(data)
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)
  if(control$verbose) cat(" Model: Additive PM-CDM \n")
  control$sampler = match.arg(control$sampler, c("ULA","MALA","RWMH"))
  control$algorithm = match.arg(control$algorithm, c("GD","ADAM","mixed"))

  if(!is.null(start.par) && !is.list(start.par) && (start.par == "random")){
    if(control$verbose) cat(" Generating random starting values for model parameters ...")
    control$iden.R = T
    if(!is.null(control$seed)) set.seed(control$seed)
    pp = pr_param_aCDM(p,q,T,control)
    if(control$verbose) cat("\r Generating random starting values for model parameters ... (Done!) \n")
  }
  if(is.null(start.par) || is.list(start.par)){
    pp = pr_param_aCDM(p,q,F,control)
    if(is.list(start.par)){
      if(!all(names(start.par) %in% names(pp))) stop("Argument `start.par' must be a list with elements `G', `mu', and/or `R'.")
      for(ii in names(start.par)){
        if(ii == "G"){
          if(nrow(start.par[[ii]]) != ncol(data)) stop("Matrix `start.par$G' mismatch rows with p.")
          start.par[[ii]] = start.par[[ii]] * cbind(1,control$Qmatrix)
        }
        if(ii == "mu") if(length(start.par[[ii]]) != q) stop("Lenght of `start.par$mu' mismatch with K.")
        if(ii == "R"){
          if((nrow(start.par[[ii]]) != ncol(start.par[[ii]])) || (max(nrow(start.par[[ii]]), ncol(start.par[[ii]])) != q)){
            stop("Rows or cols of `start.par$R' mismatch with K.")
          }
        }
        pp[[ii]] = start.par[[ii]]
      }
    }
  }

  if(!is.null(control$start.zn) && is.matrix(control$start.zn)){
    zn = control$start.zn
  } else if(control$start.zn == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn = mvtnorm::rmvnorm(nrow(data),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressWarnings(psych::fa(r = data, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                     missing = T))
    zn = tmp$scores
    if(control$cor.R) pp$R = cor(zn) else pp$R = cov(zn)
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  if(!is.null(control$seed)) set.seed(control$seed)
  out <- apmCDM_fit_rcpp(Y = data[], G = pp$G[], Qmatrix = control$Qmatrix[], mu = pp$mu[], R = pp$R[], Z= zn[], control = control)
  colnames(out$PI) <- rownames(out$G) <- colnames(data)
  colnames(out$G) <- c("(Intercept)",paste0("Z",1:q))
  colnames(out$Z) <- colnames(out$R) <- paste0("Z",1:q)
  colnames(out$U) <- paste0("U",1:q)
  if(control$return.trace){
    colnames(out$cdllk.trace) <- c("fz","fyz")
    Gnames <- paste0("G",apply(expand.grid(paste0("j",1:p),paste0("k",0:q)),1,paste,collapse = "."))
    Gnames <- Gnames[which(cbind(1,control$Qmatrix) != 0)]
    Mnames <- paste0("mu",1:q)
    Rnames <- paste0("R",apply(which(lower.tri(diag(q),diag = !control$cor.R) == T,arr.ind = T),1,paste0,collapse = ""))
    colnames(out$theta.trace) <- c(Gnames,Mnames,Rnames)
  }
  class(out) = c("apmCDM", "pmCDM")
  return(out)
}


#' Simulate data from an Partial Mastery Additive Cognitive Diagnosis Model (PM-ACDM).
#'
#' @param n Number of simulated entries.
#' @param p Number of observed (binary) variables.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param start.par (optional) List of size 3 with starting model parameters (G, mu, R).
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{Y}: A matrix \code{(n x p)} of simulated observed (binary) variables.
#'  \item \code{Z}: A matrix \code{(n x q)} of simulated latent variables (on the continuous scale).
#'  \item \code{U}: A matrix \code{(n x q)} of simulated latent variables (on the \code{[0,1]} scale).
#'  \item \code{PI}: A matrix \code{(n x p)} of predicted probabilities.
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
apmCDM_sim <- function(n, p, q, control = list(),
                        start.par = NULL, ...){
  control = pr_controlsim_aCDM(control,q,...)
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)
  if(!is.null(start.par) & !is.list(start.par) & (length(start.par) != 3))
    stop("Argument `start.par' needs to be a list with elements `G', `mu', and `R'.")
  if(!is.null(start.par)){
    pp = start.par
    if(nrow(pp$G) != nrow(control$Qmatrix)) stop("Matrix `start.par$G' mismatch rows with `control$Qmatrix'.")
    pp$G = pp$G * cbind(1,control$Qmatrix)
    if(length(pp$mu) != nrow(pp$R)) stop("Lenght of `start.par$mu' and nrows of `start.par$R' differ.")
  } else {
    pp = pr_param_aCDM(p,q,T,control)
  }
  if(!is.null(control$seed)) set.seed(control$seed)
  out <- apmCDM_sim_rcpp(n,pp$G,pp$mu,pp$R)
  out$G <- pp$G
  out$mu <- pp$mu
  out$R <- pp$R
  out$posMu <- matrix(colMeans(out$Z),nrow = n, ncol = q, byrow = T)
  out$posR <- array(cov(out$Z),dim = c(q,q,n))
  rownames(out$G) <- colnames(out$Y) <- colnames(out$PI) <- paste0("Y",1:p)
  colnames(out$G) <- c("(Intercept)",paste0("Z",1:q))
  colnames(out$R) <- rownames(out$R) <- colnames(out$posR) <- rownames(out$posR) <- colnames(out$Z) <- names(out$mu) <- colnames(out$posMu) <-paste0("Z",1:q)
  colnames(out$U) <- paste0("U",1:q)
  class(out) = c("apmCDM", "pmCDM")
  return(out)
}

#' Find number of latent variables (q) using cross-validation (PM-ACDM)
#'
#' @param Ytrain Matrix with observed binary entries to train the model. Missing entries (\code{Ytest}) coded as \code{NA}.
#' @param Ytest Matrix with observed binary entries to test the model. All entries are \code{NA} but the missing in \code{Ytrain}.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{CV.error}: Cross-validation error.
#'  \item \code{AUC}: Area under the curve error.
#'  \item \code{llk.gapmCDM} Marginal log-likelihood for the gapmCDM model
#' }
#' @details Define CV.error.
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
apmCDM_findqCV <- function(Ytrain, Ytest, q, control = list(), start.par = NULL, ...){

  control = pr_control_aCDM(control,...)
  if(control$verbose) cat(" Model: Additive PM-CDM \n")
  control$sampler = match.arg(control$sampler, c("ULA","MALA","RWMH"))
  p = ncol(data)
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)

  if(is.list(start.par)){
    if(length(start.par) != 3) stop("Argument `start.par' needs to be a list with elements `G', `mu', and `R'.")
    ppD1 = start.par
    if(nrow(ppD1$G) != ncol(data)) stop("Matrix `start.par$G' mismatch rows with p.")
    if(length(ppD1$mu) != nrow(ppD1$R)) stop("Lenght of `start.par$mu' and nrows of `start.par$R' differ.")
  }
  if(!is.null(start.par) && (start.par == "random")){
    if(control$verbose) cat(" Generating random starting values for model parameters ...")
    control$iden.R = T
    if(!is.null(control$seed)) set.seed(control$seed)
    ppD1 = pr_param_aCDM(p,q,T,control)
    if(control$verbose) cat("\r Generating random starting values for model parameters ... (Done!) \n")
  }
  if(is.null(start.par)) {
    ppD1 = pr_param_aCDM(p,q,F,control)
  }

  if(!is.null(control$start.zn) && is.matrix(control$start.zn)){
    zn = control$start.zn
  } else if(control$start.zn == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn = mvtnorm::rmvnorm(nrow(data),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressWarnings(psych::fa(r = data, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                     missing = T))
    zn = tmp$scores
    if(control$cor.R) ppD1$R = cor(zn) else ppD1$R = cov(zn)
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  if(!is.null(control$seed)) set.seed(control$seed)
  if(control$verbose) cat(paste0("\n [Training data] \n"))
  fit1 = apmCDM_fit_rcpp(Y = Ytrain[], G = ppD1$G[], Qmatrix = control$Qmatrix[], mu = ppD1$mu[], R = ppD1$R[], Z= zn[], control = control)
  out1 = apmCDM_cv_rcpp(Ytrain[], Ytest[], G = fit1$G[], mu = fit1$mu[],R = fit1$R[], Z = fit1$Z[], control = control)
  pred = ROCR::prediction(out1$Yhat,out1$Yobs)
  aucCV = unlist(methods::slot(ROCR::performance(pred, "auc"), "y.values"))

  return(list(CV.error = out1$CV.error, AUC = aucCV, llk.train = fit1$llk,
              mod.train = fit1))
}


#' Find number of latent variables (q) using cross-validated MLLK (PM-ACDM)
#'
#' @param Ytrain Matrix with observed binary entries to train the model. Missing entries (\code{Ytest}) coded as \code{NA}.
#' @param Ytest Matrix with observed binary entries to test the model. All entries are \code{NA} but the missing in \code{Ytrain}.
#' @param q Number of latent variables.
#' @param control List of control parameters (see 'Details').
#' @param ... Further arguments to be passed to \code{control}.
#'
#' @return A list with components:
#' \itemize{
#'  \item \code{mllk.test}: Test data marginal log-likelihood.
#'  \item \code{mllk.train}: Train data marginal log-likelihood.
#'  \item \code{train.mod} Fitted model on Train data (for traceability).
#' }
#' @details Define CV.error.
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
apmCDM_mllkCV <- function(Ytrain, Ytest, q, control = list(), ...){
  p = ncol(Ytrain)
  control = pr_control_aCDM(control,...)
  if(control$verbose) cat(" Model: Additive PM-CDM \n")
  if(control$verbose) cat(paste0("\n [ Training data ] \n\n"))
  control$sampler = match.arg(control$sampler, c("ULA","MALA","RWMH"))
  if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)

  ppD1 = pr_param_aCDM(p,q,F,control)

  if(!is.null(control$start.zn) && is.matrix(control$start.zn)){
    zn = control$start.zn
  } else if(control$start.zn == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn = mvtnorm::rmvnorm(nrow(Ytrain),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressMessages(suppressWarnings(psych::fa(r = Ytrain, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                                      missing = T)))
    zn = tmp$scores
    if(control$cor.R) ppD1$R = cor(zn) else ppD1$R = cov(zn)
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  Apat = as.matrix(expand.grid(lapply(1:q,function(x) c(0,1))))
  if(!is.null(control$seed)) set.seed(control$seed)
  fit1 = apmCDM_fit_rcpp(Y = Ytrain[], G = ppD1$G[], Qmatrix = control$Qmatrix[], mu = ppD1$mu[], R = ppD1$R[], Z= zn[], control = control)
  if(control$return.trace){
    colnames(fit1$cdllk.trace) <- c("fz","fyz")
    Gnames <- paste0("G",apply(expand.grid(paste0("j",1:p),paste0("k",0:q)),1,paste,collapse = "."))
    Gnames <- Gnames[which(cbind(1,control$Qmatrix) != 0)]
    Mnames <- paste0("mu",1:q)
    Rnames <- paste0("R",apply(which(lower.tri(diag(q),diag = !control$cor.R) == T,arr.ind = T),1,paste0,collapse = ""))
    colnames(fit1$theta.trace) <- c(Gnames,Mnames,Rnames)
  }

  if(control$verbose) cat(paste0("\n\n [ Testing data ] \n\n"))

  if(!is.null(control$start.zn.test) && is.matrix(control$start.zn.test)){
    zn.test = control$start.zn.test
  } else if(control$start.zn.test == "random") {
    if(control$verbose) cat(" Generating random starting values for latent variables ...")
    zn.test = mvtnorm::rmvnorm(nrow(Ytest),mean = rep(0,q))
    if(control$verbose) cat("\r Generating random starting values for latent variables ... (Done!) \n")
  } else if(control$start.zn.test == "fa") {
    if(control$verbose) cat(" Generating starting values for latent variables via Factor Analysis ...")
    tmp = suppressMessages(suppressWarnings(psych::fa(r = Ytest, nfactors = q, cor = "tet", fm = "ml", rotate = "oblimin",
                                                      missing = T)))
    zn.test = tmp$scores
    if(control$verbose) cat("\r Generating starting values for latent variables via Factor Analysis ... (Done!) \n")
  }

  if(!control$cv.useFitPos){
    IS_posMu = matrix(c(fit1$mu), nrow = nrow(Ytest), ncol = q)
    IS_posR = array(fit1$R, dim = c(q,q,nrow(Ytest)))
  } else {
    IS_posMu = fit1$posMu
    IS_posR =  fit1$posR
  }

  testmllk = fy_aCDM_IS(Ytest[],G = fit1$G[],
                        mu = fit1$mu[],L = t(chol(fit1$R[])), Z = zn.test[],
                        pM = IS_posMu[], pR = IS_posR[], control = control)
  return(list(mllk.test = testmllk, mllk.train = fit1$llk, train.mod = fit1)) #
}
