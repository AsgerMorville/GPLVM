#Generate data function.
library(MASS)
generateData <- function(t,p,q,parameters){
  tvec <- seq(from=1,to=t,length.out = t)
  Kt <- fillup(as.matrix(tvec),matern,list(l=20,r=1))
  x<-mvrnorm(n=1,mu=rep(0,t),Sigma=Kt)
  Ktinv <- solve(Kt)
  Kx <- fillup(as.matrix(x),ardkernel,list(sigmaf=1,lengthscales=1))
  y <- matrix(NA,nrow=p,ncol=t)
  f <- y
  for (i in 1:p){
    f[i,] <- mvrnorm(n=1,mu=rep(0,t),Sigma=Kx)
    y[i,] <- rnorm(t,mean=(f[i,]),sd=0.1)
  }
  y <- t(y)
  f <- t(f)
  return(list(y=y,f=f,x=x))
}

