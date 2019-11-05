#Model to contain S,  fhatk,Ak, t, Kx, 
library(R.utils)
library(MASS)
sourceDirectory("Modules/",modifiedOnly = F)
#simulate data
model <- list()

model$t <- 30
model$q <- 2
model$p <- 20
model$tmat <- matrix(1:model$t,ncol=1)
model$thetax <- list(r=1,l=1)
model$thetaf <- list(rho=1,delta=1)
Kt <- fillup(model$tmat,matern,model$thetax)
#Draw an x-sample from this
model$Ktinv <- solve(Kt)
x1<- mvrnorm(n=1,mu=rep(0,model$t),Sigma = Kt)
x2<- mvrnorm(n=1,mu=rep(0,model$t),Sigma = Kt)
model$x <- matrix(c(x1,x2),ncol=model$q,nrow=model$t)

#draw f and y-samples
Kx <- fillup(model$x,ARDKernel,model$thetaf)

y <- matrix(NA,nrow=model$t,ncol=model$p)
f <- matrix(NA,nrow=model$t,ncol=model$p)
for (i in 1:model$p){
  f[,i] <- mvrnorm(n=1,mu=rep(0,model$t),Sigma=Kx)
  y[,i] <- rpois(model$t,exp(f[,i]))
}
model$y <-y
model$fhatk <- f



#Prepare for objfunc
Kx_object <- Kx_obj(model$x)
A_arr <- array(NA,dim=c(model$t,model$t,model$p))
for (i in 1:model$p){
  A_arr[,,i] <- diag(exp(f[,i]))-Kx_object$Kxinv
}
model$A_arr <- A_arr
Sk_arr <- array(NA,dim=c(model$t,model$t,model$p))
Skinv_arr <- array(NA,dim=c(model$t,model$t,model$p))
for (i in 1:model$p){
  Skinv_arr[,,i] <- A_arr[,,i]-Kx_object$Kxinv
  Sk_arr[,,i] <- solve(Skinv_arr[,,i])
}
model$Skinv_arr <- Skinv_arr
model$Sk_arr <- Sk_arr
#model$Kx <- Kx_obj$Kx
#model$Kx_in

#Can we optimse?

start <- eigenInit(y,2)

plot(start[,1])
plot(model$x[,1])
test <- scg(objfunc,gradObj,as.vector(start),1)


#Testing
source("Test/objfuncc.R")
objfunc(as.vector(model$x))
library(numDeriv)
test <- rnorm(60)
num <- grad(objfunc,as.vector(model$x))
reel<-gradObj(as.vector(model$x))
reel
plot(num-as.vector(reel))

num2 <- grad(objfunc,as.vector(test))
reel2 <- gradObj(as.vector(test))
plot(num2-reel2)

num
plot(as.vector(reel))
plot(num)
