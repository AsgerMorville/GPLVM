a <- as.matrix(read.table("skildpaddedata"))
#Main script
library(MASS)
library(R.utils)
sourceDirectory("GaussianReg/GPLVM-trad/Modules/")
set.seed(1)
theta <- c(1,1,0.1)
n <- 45
t <- 39
q <- 1
yDim <- n
tvec <- seq(from=1,to=39,length.out = 39)
matern <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}
Kt <- fillup(tvec,tvec,matern,1)
#x <- matrix(NA,nrow=1,ncol=t)
x<-mvrnorm(n=1,mu=rep(0,t),Sigma=Kt)
x <- scale(x,scale=F)
Ktinv <- solve(Kt)
plot(x)
l <- 1
par <- c(1,1,0.1)
sigma <- 0.1
#x <- seq(from=0.001,to=1,length.out=n)
Kx <- fillup2(as.matrix(x),RBF_vec,par)
#diag(t)*sigma^2

y <- matrix(NA,nrow=n,ncol=t)
f <- y
for (i in 1:n){
  f[i,] <- mvrnorm(n=1,mu=rep(0,t),Sigma=Kx)
  y[i,] <- rnorm(t,mean=(f[i,]),sd=0.1)
}


y <- t(y)
f <- t(f)

plot(x,y[,2])
#Data generated


testfunc <- function(x){
  return(logLik(c(x[1],x[2],x[3]),as.matrix(x[4:42])))
}
testgrad <- function(x){
  X <- as.matrix(x[4:42])
  obj <- gradLogLik(c(x[1],x[2],x[3]),X)
  hyper <-obj$hyper
  x <- as.vector(obj$x)
  return(c(hyper,x))
}
#testgrad(c(1,1,0.1,x))

#logLik(c(1,1,0.1),as.matrix(x))

gradLogLik(c(1,1,0.1),as.matrix(x))

test <- scg(testfunc,testgrad,c(1,1,0.1,x),1)
test2 <- scg(testfunc,testgrad,c(2,2,0.2,100*eigenInit(y,1)),1)
plot(1:39,test2[4:42])


plot(x,y[,2])
plot(test2[4:103],y[,2])

plot(tvec,x)
plot(tvec,test2[4:103])
plot(tvec,eigenInit(y,1))
