#Main script
library(MASS)
library(R.utils)
sourceDirectory("GaussianReg/GPLVM-trad/Modules/",modifiedOnly = F)
set.seed(1)
theta <- c(1,1,0.1)
n <- 40
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

plot(x,y[,7])
plot(x[,1],x[,2])
#Data generated


### Crete model object
modeltest <- list()

modeltest$x <- x
modeltest$y <- y
modeltest$hyper$gamma <- 1
modeltest$hyper$alpha <- 1
modeltest$YYt <- y%*%t(y)
modeltest$Ktinv <- Ktinv
modeltest$hyper$beta <- 1

###
vecGrad(c(1,1,100,x))


logLikv2(modeltest)
gradLogLikv2(modeltest)
####

startx <- eigenInit(y,1)
plot(x)
plot(startx)
test <- scg(vecLogLik,vecGrad,c(1,1,100,startx),1)



plot(x)
plot(test[-c(1:3)])

#Test gradients
testfunc <- function(vector){
  modeltest$hyper$gamma <- vector[1]
  modeltest$hyper$alpha <- vector[2]
  modeltest$hyper$beta <- vector[3]
  modeltest$x <- as.matrix(vector[-c(1:3)])
  return(logLik(modeltest))
}

testfunc(c(1,1,10,x))

ourguess <- function(vector){
  modeltest$hyper$gamma <- vector[1]
  modeltest$hyper$alpha <- vector[2]
  modeltest$hyper$beta <- vector[3]
  modeltest$x <- as.matrix(vector[-c(1:3)])
  obj <- gradLogLikv2(modeltest)
  
  return(c(obj$hyper,obj$x))
}

num <- grad(testfunc,c(1,1,10,x))
reel <- ourguess(c(1,1,10,x))
plot(num-reel)
ourguess(c(1,1,10,x))

ourguess <- function(betaa){
  modeltest$hyper$beta <- betaa
  return(gradLogLikv2(modeltest)$hyper[3])
}
ourguess(10)

grad(testfunc,10)
###

testfunc(100)

x <- model$x
y <- model$y
D <- dim(y)[2]
gamma <- model$hyper$gamma
alpha <- model$hyper$alpha
Ktinv <- model$Ktinv
YYt <- model$YYt
##



testfunc <- function(x){
  return(logLik(c(x[1],x[2],x[3]),as.matrix(x[4:103])))
}
testgrad <- function(x){
  X <- as.matrix(x[4:103])
  obj <- gradLogLik(c(x[1],x[2],x[3]),X)
  hyper <-obj$hyper
  x <- as.vector(obj$x)
  return(c(hyper,x))
}
#testgrad(c(1,1,0.1,x))

#logLik(c(1,1,0.1),as.matrix(x))

gradLogLik(c(1,1,0.1),as.matrix(x))

test <- scg(testfunc,testgrad,c(1,1,0.1,x),1)
test2 <- scg(testfunc,testgrad,c(2,2,0.2,eigenInit(y,1)),1)

plot(x,y[,2])
plot(test2[4:103],y[,2])

plot(tvec,x)
plot(tvec,test2[4:103])
plot(tvec,eigenInit(y,1))
