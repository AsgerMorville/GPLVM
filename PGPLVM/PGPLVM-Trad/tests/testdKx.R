#Tests

library(numDeriv)

num <- grad(testfunction,Xny)
reel <- vecGrad(Xny)

plot(num)
#Test of dKxdX

ourguess <- function(s){
  x[3,1] <- s
  Kx <- fillup(x,x,kern,1)
  mat <- difkernq(3,1,x,Kx,1)
  return(as.vector(mat))
}

matrix(ourguess(3),ncol=100,nrow=100)

#The numderiv
numerr <- function(s){
  x[3,1] <- s
  Kx <- fillup(x,x,kern,1)
  return(as.vector(Kx))
}
numer <- jacobian(numerr,1)
reel <- ourguess(1)
matrix(reel-numer,ncol=100,nrow=100)
plot(numer-reel)
