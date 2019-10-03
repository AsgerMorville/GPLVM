gradLoglik <- function(model){
  #precomputations, obtain object mod containing
  #tvec, y, mu, S, xu,hyperparams and
  #psistatistics, Kuu, Kuuinv, Kxx, Kxxinv
  
  mod <- modCalculator(model)
  #With the mod object, we can compute derivatives
  #wrt mu, S, xu and hyper
  #dFdMubar-matrix
  
  gradobj <- list()
  
  gradobj$mubar <- dFdmubar(mod)
  gradobj$Lambda  <- dFdlambda(mod)
  #gradobj$xu <- dFdxu(mod)
  #gradobj$thetaf <- dFdthetaf(mod)
  #gradobj$thetax <- dFdthetax(mod)
  #gradobj$beta <- dFdbeta(mod)
  return(gradobj)
}
