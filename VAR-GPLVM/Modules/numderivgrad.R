#numderiv
#?grad
grad(vecLogLik,startguess)
vecGradLogLik(startguess)
numderivgrad <- function(x){
  return(grad(vecLogLik,x))
}