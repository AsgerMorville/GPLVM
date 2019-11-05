matern2 <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}