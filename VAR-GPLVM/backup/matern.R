matern <- function(dif,par){
  return(par$r*exp(-dif/par$l))
  #return(theta[2]^2*exp(-r^2/(2*theta[1]^2)))
}
