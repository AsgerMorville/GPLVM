ardkernel <- function(dif,theta){
  return(theta$sigmaf*exp(-0.5*sum(dif^2/theta$lengthscales)))
  #return(theta[2]^2*exp(-r^2/(2*theta[1]^2)))
}
