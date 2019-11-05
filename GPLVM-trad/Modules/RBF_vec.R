#RBF kernel
RBF_vec <- function(r,theta){
  return(theta$alpha*exp(-theta$gamma*r/2))
}
