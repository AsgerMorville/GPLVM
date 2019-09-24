#Kernel derivatives

dKl <- function(r,theta){
  return(theta[2]^2*exp(-r/(2*theta[1]^2))*(-r/(theta[1]^3)))
}
dKa <- function(r,theta){
  return(2*theta[2]*exp(-r/(2*theta[1]^2)))
}
dKsigma <- function(sigma){
  return(diag(n)*sigma*2)
}

