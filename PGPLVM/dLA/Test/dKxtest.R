##Tests

testfunc <- function(x3){
  x <- model$x
  x[27,2] <- x3
  Kx <- fillup(x,ARDKernel,theta=model$thetaf)
  return(as.vector(Kx))
}
ourgrad <- function(x3){
  x <- model$x
  x[27,2] <- x3
  Kx <- fillup(x,ARDKernel,theta=model$thetaf)
  #print(Kxtemp$Kx)
  return(dKxdxij(2,2,Kx,model$thetaf,x))
}

#library(numDeriv)
num <- jacobian(testfunc,0.1)
reel <- as.vector(ourgrad(0.1))

plot(reel)

#model$thetaf$delta
#plot(num)