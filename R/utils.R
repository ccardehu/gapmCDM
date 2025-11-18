#' Compute marginal log-likelihood (via Importance Sampling) for an object of class "pmCDM".
#'
#' @param mod \code{pmCDM} object.
#' @param Y Dataset to evaluate MLLK (default \code{NULL})
#' @param control List for importance sampling options (\code{nsim})
#'
#' @return A double:
#' \itemize{
#'  \item \code{llk}: Marginal log-likelihood (computed via Monte-Carlo integration) on \code{Y} evaluated at estimated parameters.
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
logLik.pmCDM <- function(mod, Y = NULL, control = list(), ...){
  if(!is.null(mod$Y)) Y <- as.matrix(mod$Y)
  if(!is.null(Y) && !is.matrix(Y)) Y <- as.matrix(Y)
  if("gapmCDM" %in% class(mod)){
    if(!is.null(mod$llk)){
      return(mod$llk)
    } else {
      control = pr_control_gaCDM(control,...)
      control$basis = match.arg(control$basis, c("is","bs","pwl"))
      if(is.null(control$degree) && control$basis == "pwl") control$degree = 1
      if(is.null(control$degree) && control$basis != "pwl") control$degree = 2
      if(control$verbose) cat(" Model: Generalized Additive PM-CDM")
      llk = fy_gapmCDM_IS(Y = Y[],A = mod$A[],C = mod$C[],mu = mod$mu[],L = t(chol(mod$R[])),Z = mod$Z[],
                          pM = mod$posMu[], pR = mod$posR[], control = control)
      return(llk)
    }
  }
  else if("apmCDM" %in% class(mod)){
    if(!is.null(mod$llk)){
      return(mod$llk)
    } else {
      control = pr_control_aCDM(control,...)
      if(control$verbose) cat(" Model: Additive PM-CDM")
      q = ncol(mod$R)
      p = ncol(Y)
      if(is.null(control$Qmatrix)) control$Qmatrix = matrix(1,p,q)
      llk = fy_aCDM_IS(Y = Y[],G = mod$G[],
                       mu = mod$mu[],L = t(chol(mod$R[])),Z = mod$Z[],
                       pM = mod$posMu[], pR = mod$posR[], control = control)
      return(llk)
    }
  }
}


#' Compute fitted probabilities for an object of class "pmCDM".
#'
#' @param mod \code{pmCDM} object.
#' @param U Points for latent attributes on which compute the fitted probabilities
#' @param control Additional values for control
#'
#' @return A matrix with fitted probabilities (one entry per row in \code{U}):
#' \itemize{
#'  \item \code{PI}: Matrix of fitted probabilities.
#' }
#' @details Test
#' @author Camilo Cárdenas-Hurtado (\email{c.a.cardenas-hurtado@@lse.ac.uk}).
#' @export
fitted_probs <- function(mod, U, control = list(), ...){
  if(!is.null(U) && !is.matrix(U)) U <- as.matrix(U)
  if("gapmCDM" %in% class(mod)){
    control = pr_control_gaCDM(control,...)
    control$basis = match.arg(control$basis, c("is","bs","pwl"))
    if(is.null(control$degree) && control$basis == "pwl") control$degree = 1
    if(is.null(control$degree) && control$basis != "pwl") control$degree = 2
    if(control$verbose) cat(" Model: Generalized Additive PM-CDM")
    if(control$basis != "is"){
      SpUObj = SpU_bsp(U,control$knots,control$degree)
      isM = SpUObj[,,1]
    } else {
      SpUObj = SpU_isp(U,control$knots,control$degree)
      isM = SpUObj[,,1]
    }
    if(ncol(mod$C) != (length(control$knots) + 1)) stop("Check input knots (add argument control$knots congruent with mod)")
    PI = prob(mod$A[],mod$C[],isM)
    return(PI)
  }
  else if("apmCDM" %in% class(mod)){
    control = pr_control_aCDM(control,...)
    if(control$verbose) cat(" Model: Additive PM-CDM")
    q = ncol(mod$R)
    PI = prob_aCDM(mod$G[],U)
    return(PI)
  }
}
