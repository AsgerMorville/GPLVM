#RBF kernel
RBF_vec <- function(r,theta){
  return(theta[2]^2*exp(-r/(2*theta[1]^2)))
}
