gradientNewton <- function(fvec,y,Kx_inv){
  return(as.vector(y-exp(fvec)-Kx_inv%*%fvec))
}
