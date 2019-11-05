matern <- function(r,thetax){
  return(thetax$r*exp(-sqrt(r)/(thetax$l)))
}
