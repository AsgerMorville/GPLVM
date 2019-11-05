#ARDKernel
ARDKernel <- function(r,theta){
  return(theta$rho*exp(-r/(2*theta$delta^2)))
}

