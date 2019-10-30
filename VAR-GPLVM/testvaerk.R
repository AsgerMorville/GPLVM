#numderiv
#?grad
library(numDeriv)

num<-grad(vecLogLik,end)
reel <- vecGradLogLik(end)
num<-grad(vecLogLik,startguess)
reel <- vecGradLogLik(startguess)
numderivgrad <- function(x){
  return(grad(vecLogLik,x))
}

fillup(matrix(1:10),ardkernel,thetaf)

muu = matrix(1:n)
muu[2] <- 10
(matrix(muu[,1],n,n,byrow=F)-matrix(muu[,1],n,n,byrow=T))
plot(num[1:60])
plot(reel[1:60])
plot(num[61:90])
plot(reel[61:90])
plot(num[31:60]-reel[31:60])
plot(num[1:59]-reel[1:59])
vecLogLik(startguess)
plot(num)
plot(reel)
plot(num-reel)
startguess[91] <- 1
