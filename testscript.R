library(MASS)
#Create example
#F is nonlinear, eg. sin
q <- 1
dim <-20
n <- 100
x <- matrix(runif(q*n,min=-5,max=5),nrow=n,ncol=q)
#x<- rnorm(n,mean=0,sd=2)
l <- 1
sigma <- 0.1
#x <- seq(from=0.001,to=1,length.out=n)
Ky <- fillup2(x,kern,c(l,1,0.1))+diag(n)*sigma^2

y <- matrix(NA,nrow=n,ncol=dim)
for (i in 1:dim){
  y[,i] <- mvrnorm(n=1,mu=rep(0,n),Sigma=Ky)
}
y <- scale(y,scale=F)
plot(y[,1],y[,12])


### unwrapper
testfunc <- function(x){
  X <-matrix(x,ncol=q,nrow=n,byrow=F)
  return(logLik(c(1,1,0.1),X))
}
testgrad <- function(x){
  X <-matrix(x,ncol=q,nrow=n,byrow=F)
  obj <- gradLogLik(c(1,1,0.1),X)
  return(as.vector(obj))
}

test2 <- scg(testfunc,testgrad,x,0.1)

test3 <- scg(testfunc,testgrad,Xny,0.5)
test3 <- scg(testfunc,testgrad,runif(100),0.1)
testgrad(runif(100))
###
test
plot(y[,1],y[,2])
plot(x,y[,1],xlim=c(-4,4))
plot(test2,y[,1],xlim=c(-4,4))
plot(x,y[,2])
testfunc(x)
testfunc(test3)
testfunc(runif(100))
y <- scale(y,scale=F)
Test <- scg(loglikchol,gradientChol2.0,c(1,1,0.1,Xny),1)

YYt <- y%*%t(y)

alpha <- 1
beta <- 1

eig <- eigen(YYt)
D <- dim
lj <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
Xny <- eig$vectors[,1]*lj

plot(Xny,y[,2])
plot(x,y[,2],xlim = c(-4,4))

test2 <- optim(Xny,loglikchol,gradientchol,method="CG")
test3 <- scg(loglikNeil,gradNeil,Xny,0.1)
plot(test$par,y[,2],xlim=c(-4,4))
l <- 1
sigma <- 0.1
test <- optim(Xny,loglikchol,gradientchol,method="CG")



l <- 1
sigma <- 0.1
loglikchol(x)
xpre <- runif(100,min=-2.5,max=2.5)
realy <- predictor(xpre,x)
fakey <- predictor(xpre,test$par)
realfake <- predictor(xpre,Xny)
plot(-xpre,fakey)
plot(xpre,realy)
plot(-xpre,realfake)

realy2 <- predictor(xpre,x,2)
fakey2 <- predictor(xpre,test$par,2)
realfake2 <- predictor(xpre,Xny,2)
plot(xpre,realy2)
plot(-xpre,fakey2)
plot(-xpre,realfake2)
Ysc <- scale(Y,scale = F)
plot(Ysc[,1],Ysc[,2])
YYt <- Ysc%*%t(Ysc)
alpha <- 1
beta <- 1

gradientchol(-xpre)

eig <- eigen(YYt)
D <- 2
lj <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
Xny <- eig$vectors[,1]*lj

plot(Xny)
plot(x)
l <-1
sigma <- 0.01
optim(Xny,loglikchol,gradientchol,method="BFGS")
gradientchol(Xny)


Test<-optim(Xny,loglikchol,gradientchol,method="BFGS")
Test2 <- optim(Xny,loglikchol,gradientchol,method="CG")
Test3 <- optim(Xny,loglikchol,gradientchol,method="CG", control = list(type = 3))
Test22 <- optim(Xny,loglikchol,gradientchol,method="CG", control = list(type = 2))
