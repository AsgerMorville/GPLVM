#Approximated Laplace 
set.seed(3)
library(MASS)
source("PGPLVM/PGPLVM-Trad/Modules/funks.R")
y <- as.matrix(read.table("skildpaddedata"))
#theta <- c(1,1,0.1)
n <- 45
t <- 39
q <- 1
yDim <- n
tvec <- 1:t
matern <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}
Kt <- fillup(tvec,tvec,matern,20)
#x <- matrix(NA,nrow=1,ncol=t)
x<-mvrnorm(n=1,mu=rep(0,t),Sigma=Kt)
Ktinv <- solve(Kt)
plot(x)

x <- scale(x,scale=F)
l <- 1
par <- c(1,1,0.1)
sigma <- 0.1
#x <- seq(from=0.001,to=1,length.out=n)
Kx <- fillup(x,x,kern,l)
#diag(t)*sigma^2

y <- matrix(NA,nrow=n,ncol=t)
f <- y
for (i in 1:n){
  f[i,] <- mvrnorm(n=1,mu=rep(0,t),Sigma=Kx)
  y[i,] <- rpois(t,exp(f[i,]))
}
plot(tvec,y[11,])
lines(tvec,exp(f[11,]))
aLA(y,x,250)
plot(y[4,],y[5,])
plot(x,f[4,])

aLA(tvec,y,x,250)
?rpois


testfunction <- function(x){
  xmat <- as.matrix(x)
  fullLogLik(tvec,y,xmat,250)
}
testgrad <- function(x){
  xmat <- as.matrix(x)
  return(as.vector(fullGrad(tvec,y,xmat,250)))
}

y <- t(y)
alpha <- 1
beta <- 1

eig <- eigen(t(y)%*%(y))
D <- 2
lj <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
Xny <- eig$vectors[,1]*lj

plot(1:39,Xny)
test<-scg(testfunction,testgrad,100*Xny,1)

plot(test,)
plot(Xny,x)
plot(tvec,test)
plot(tvec,x)
plot(tvec,test)
plot(x,f[2,])
plot(-realx2,f[2,])
plot(-realx2,f[2,],xlim=c(-3,3))
plot(x,f[2,],xlim = c(-3,3))
