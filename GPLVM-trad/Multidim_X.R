#Multidim X
#Main script
library(MASS)
library(R.utils)
sourceDirectory("GaussianReg/GPLVM-trad/Modules/",modifiedOnly=FALSE)

thetaf <- list(alpha=1,gamma=1)
n <- 40

t <- 100
q <- 2
yDim <- n
tvec <- seq(from=1,to=t,length.out = t)
matern <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}
Kt <- fillup(tvec,tvec,matern,20)
#x <- matrix(NA,nrow=1,ncol=t)
x <- matrix(NA,t,2)
x[,1]<-mvrnorm(n=1,mu=rep(0,t),Sigma=Kt)
x[,2]<-mvrnorm(n=1,mu=rep(0,t),Sigma=Kt)
x <- scale(x,scale=F)
Ktinv <- solve(Kt)
plot(x[,1])
plot(x[,2])
plot(x[,1],x[,2])
#l <- 1
#par <- c(1,1,0.1)
#sigma <- 0.1
#x <- seq(from=0.001,to=1,length.out=n)
Kx <- fillup2(as.matrix(x),RBF_vec,list(alpha=1,gamma=1))
#Kx
x#diag(t)*sigma^2
#Kx
y <- matrix(NA,nrow=n,ncol=t)
f <- matrix(NA,nrow=n,ncol=t)
for (i in 1:n){
  f[i,] <- mvrnorm(n=1,mu=rep(0,t),Sigma=Kx)
  y[i,] <- rnorm(t,mean=(f[i,]),sd=0.1)
}


y <- t(y)
f <- t(f)

plot(x[,2],y[,4])
#Data generated
dff <- data.frame(1:100)
dff$x1 <- x[,1]
dff$x2 <- x[,2]
library(ggplot2)


qplot(x1,x2 , data = dff, colour = X1.100)

modeltest <- list()

modeltest$Xny <- eigenInit(y,2)
modeltest$y <- y
modeltest$hyper$gamma <- 1
modeltest$hyper$alpha <- 1
modeltest$YYt <- y%*%t(y)
modeltest$Ktinv <- Ktinv
modeltest$hyper$beta <- 100
modeltest$q <- dim(modeltest$Xny)[2]

#Xny <- eigenInit(y,2)
#plot(1:100,Xny[,2])
startguess <- c(1,1,100,as.vector(modeltest$Xny))
vecLogLik(startguess)

library(numDeriv)
num <- grad(vecLogLik,startguess)
reel <- vecGrad(startguess)
plot(num)
plot(-reel)
plot(reel-num)
reel
num
plot(num[1:103]-reel[1:103])
reel
test <- scg(vecLogLik,vecGrad,c(1,1,100,as.vector(modeltest$Xny)),1)
plot(Xny[,2])
plot(test[104:203])
lines(test[104:203])
plot(test[4:103])
lines(test[4:103])

plot(x[,1])
plot(-test[4:103])

vecGrad(c(1,1,100,as.vector(x)))
vecLogLik(c(1,1,100,as.vector(x)))



testfunc <- function(x){
  return(logLik(c(1,1,0.1),as.matrix(x)))
}
testgrad <- function(x){
  x <- matrix(x,ncol=q,byrow=F)
  obj <- gradLogLik(c(1,1,0.1),x)
  #hyper <-obj$hyper
  x <- as.vector(obj$x)
  return(x)
}
testgrad(x)
logLik(c(1,1,0.1),as.matrix(x))

gradLogLik(c(1,1,0.1),as.matrix(x))
test
test <- scg(testfunc,testgrad,x,1)
plot(x,y[,2])
plot(test,y[,2])

Xny <- eigenInit(y,2)
test <- scg(testfunc,testgrad,Xny,1)

plot(x[2,],y[,2],xlim=c(-1,2))
plot(tvec,test[,1])
plot(tvec,x[1,])
plot(Xny)
plot(tvec,test[,1])
plot(tvec,x[,1])
plot(x)
lines(x,y[,2])
plot(x,y[,2])
plot(tvec,x[,1])
plot(x)
plot(test)
plot(tvec,test[,1])
plot(tvec,test)
lines(tvec,x)
