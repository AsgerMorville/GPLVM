#Main script
library(MASS)
library(R.utils)
sourceDirectory("Modules/")
set.seed(1)
theta <- c(1,1,0.1)
n <- 20
t <- 100
q <- 1
yDim <- n
tvec <- seq(from=1,to=100,length.out = 100)
matern <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}
Kt <- fillup(tvec,tvec,matern,20)
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
  return(logLik(c(1,1,0.1),as.matrix(x)))
}
testgrad <- function(x){
  return(gradLogLik(c(1,1,0.1),as.matrix(x)))
}
logLik(c(1,1,0.1),as.matrix(x))

gradLogLik(c(1,1,0.1),as.matrix(x))

test <- scg(testfunc,testgrad,x,1)
plot(x,y[,2])
plot(test,y[,2])
YYt <- y%*%t(y)

alpha <- 1
beta <- 1

eig <- eigen(YYt)
D <- dim(y)[2]
lj <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
Xny <- eig$vectors[,1]*lj
test <- scg(testfunc,testgrad,Xny,1)
plot(test,y[,2],xlim=c(-1,2))
lines(x,y[,2])
plot(x,y[,2])
plot(tvec,x)
plot(tvec,test)
lines(tvec,x)
