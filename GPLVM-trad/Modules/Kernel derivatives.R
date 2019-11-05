#Kernel derivatives

dKg <- function(r,theta){
  return(theta$alpha*exp(-theta$gamma*r/2)*(-(r)/2))
}
dKa <- function(r,theta){
  return(exp(-r*theta$gamma/2))
}
dKsigma <- function(beta){
  t <- dim(modeltest$y)[1]
  return(-diag(t)/beta^2)
}

