#Testscript2
library(MASS)
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
Ktinv <- solve(Kt)
plot(x)
l <- 1
par <- c(1,1,0.1)
sigma <- 0.1
#x <- seq(from=0.001,to=1,length.out=n)
Kx <- fillup(x,x,kern1,l)
#diag(t)*sigma^2

y <- matrix(NA,nrow=n,ncol=t)
f <- y
for (i in 1:n){
  f[i,] <- mvrnorm(n=1,mu=rep(0,t),Sigma=Kx)
  y[i,] <- rnorm(t,mean=(f[i,]),sd=0.1)
}


y <- t(y)
f <- t(f)

#y <- scale(y,scale=F)
plot(tvec,y[,7])

lines((f[,7]))
plot(x)
plot(x,y[,8])
plot(y[,8],y[,9])
plot(y[3,],y[4,])

testfunc <- function(x){
  X <-matrix(x,ncol=q,nrow=t,byrow=F)
  return(logLik(X))
}
testgrad <- function(x){
  X <-matrix(x,ncol=q,nrow=t,byrow=F)
  obj <- gradLogLik(X)
  return(as.vector(obj))
}
rep(0,100)
nyx<-scg(testfunc,testgrad,x,1)
testeer <- optim(x,fn=testfunc,gr=testgrad,method="BFGS")
alpha <- 1
beta <- 1

eig <- eigen(y%*%t(y))
D <- 2
lj <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
Xny <- eig$vectors[,1]*lj

plot(Xny,y[,1])
plot(Xny)
plot(x)
estfunc(as.matrix(x))
testgrad(x)
gradLogLik(c(1,1,0.1),as.matrix(x))
nyxx <- scg(testfunc,testgrad,Xny,1)
plot(-nyxx,y[,1])
plot(Xny,y[,1])
plot(x,y[,1])
plot(x,nyxx)
