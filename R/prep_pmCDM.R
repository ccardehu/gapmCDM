pr_control_gaCDM <- function(control, ...){

  con <- list("burn.in" = 3e3, "iter.lim" = 1e4, "tune.lim" = 50,
              "stop.eps" = 1e-5,
              "gamma" = 1, "gamma.A" = 1, "gamma.CD" = 1, "gamma.R" = 1,
              "h" = 1e-2, "tune.gamma" = 1, "return.trace" = F,
              "basis" = "pwl", "degree" = NULL, "knots" = seq(.1,.9,by = 0.1),
              "Qmatrix" = NULL,
              "cor.R" = T,
              "nsim" = 1e4, "verbose" = T, "verbose.every" = 10,
              "seed" = NULL, "mu" = NULL, "R" = NULL, "sampler" = "ULA",
              "start.zn" = "fa", "start.zn.test" = "random",
              "window" = 10, "stop.atconv" = T,
              "algorithm" = "GD", "adam.b1" = .9, "adam.b2" = .999,
              "sa.power" = -.51, "cv.useFitPos" = T)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    stop("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

pr_control_aCDM <- function(control, ...){

  con <- list("burn.in" = 3e3, "iter.lim" = 1e4, "tune.lim" = 50,
              "stop.eps" = 1e-5,
              "gamma" = 1, "gamma.G" = 1, "gamma.mu" = 1, "gamma.R"= 1,
              "h" = 1e-2, "tune.gamma" = 1, "return.trace" = F,
              "max.guess" = 0.10, "max.slip" = 0.0, "Qmatrix" = NULL, "allow.slip" = T,
              "cor.R" = F,
              "nsim" = 1e4, "verbose" = T, "verbose.every" = 10,
              "seed" = NULL, "mu" = NULL, "R" = NULL, "sampler" = "ULA",
              "start.zn" = "fa", "start.zn.test" = "random",
              "window" = 10, "stop.atconv" = T,
              "damp.factor" = 1,
              "algorithm" = "GD", "adam.b1" = .9, "adam.b2" = .999,
              "sa.power" = -.51, "cv.useFitPos" = T)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    stop("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

pr_controlsim_gaCDM <- function(control,q,...){

  con <- list("degree" = NULL, "knots" = seq(.1,.9,by=0.1), "prob.sparse" = 0.75,
              "iden.R" = F, "seed" = NULL, "mu" = NULL, "R" = NULL,
              "Qmatrix" = NULL, "A.equal" = F,
              "basis" = "pw")
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    stop("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

pr_controlsim_aCDM <- function(control,q,...){

  con <- list("iden.R" = F, "max.guess" = 0.10, "max.slip" = 0.0,
              "seed" = NULL, "mu" = NULL, "R" = NULL, "Qmatrix" = NULL)
  control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(namc[!namc %in% namC]) > 0)
    stop("Unknown names in control: ", paste(namc[!namc %in% namC], collapse = ", "))
  return(con)
}

#' @export
pr_param_gaCDM <- function(p,q,tp,sim = F,control){
  if(!sim){
      As <- matrix(1/q,p,q) * control$Qmatrix
      As <- t(apply(As,1,function(x){x * 1/sum(x)}))
      # Ds <- array(1, dim = c(p,tp,q))
      Ds <- array(1/tp, dim = c(p,tp,q))
      if(control$basis == "is"){
        Cs <- array(1/tp,dim = c(p,tp,q))
      } else {
        Cs <- D2C(Ds)
      }
      Rs <- diag(q)
      mu <- rep(0,q)
      return(list("A" = As, "C" = Cs, "D" = Ds, "mu" = mu, "R" = Rs))
  } else {
    simp <- genpar(p,q,tp,control$prob.sparse,control$Qmatrix,control$basis)
    if(control$A.equal){
      As <- matrix(1/q,p,q) * control$Qmatrix
      As <- t(apply(As,1,function(x){x * 1/sum(x)}))
    } else {
      As <- simp$A
    }
    Cs <- simp$C
    Ds <- simp$D
    if(is.null(control$R)){
      if(control$iden.R){
        control$R = diag(q)
      } else {
        control$R = simstudy::genCorMat(q, cors = seq(-0.5,0.5,length.out = q*(q-1)/2))
      }
    }
    if(is.null(control$mu)){
      control$mu = rep(0,q)
    }
    return(list("A" = As, "C" = Cs, "D" = Ds, "mu" = control$mu, "R" = control$R))
  }
}

#' @export
pr_param_aCDM <- function(p,q,sim = F,control){
  if(!sim){
    Gs <- matrix(1/q,p,q) * control$Qmatrix
    Gs <- t(apply(Gs,1,function(x){x * 1/sum(x)}))
    if(q == 1){Gs <- t(Gs)}
    for(i in 1:p){
      Gs[i,Gs[i,] != 0] <- Gs[i,Gs[i,] != 0] - (control$max.guess + control$max.slip)/sum(control$Qmatrix[i,])
    }
    Gs <- cbind(control$max.guess, Gs)
    Rs <- diag(q)
    mu <- rep(0,q)
    return(list("G" = Gs, "mu" = mu, "R" = Rs))
  } else {
    Gs <- genpar_aCDM(control$Qmatrix, control$max.guess, control$max.slip)
    if(is.null(control$R)){
      if(control$iden.R){
        control$R = diag(q)
      } else {
        control$R = simstudy::genCorMat(q, cors = seq(-0.5,0.5,length.out = q*(q-1)/2))
      }
    }
    if(is.null(control$mu)){
      control$mu = rep(0,q)
    }
    return(list("G" = Gs, "mu" = control$mu, "R" = control$R))
  }
}
