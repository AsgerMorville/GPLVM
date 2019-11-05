library(MASS)
#Create example
#F is nonlinear, eg. sin
n <- 20
t <- 100

tvec <- seq(from=1,to=100,length.out = 100)
matern <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}
Kt <- fillup(tvec,tvec,matern,20)
#x <- matrix(NA,nrow=1,ncol=t)
x<-mvrnorm(n=1,mu=rep(0,t),Sigma=Kt)

plot(x)
#x <- runif(t,min=-5,max=5)
#x<- rnorm(n,mean=0,sd=2)
l <- 1
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

#y <- scale(y,scale=F)
plot(y[3,])
lines(exp(f[3,]))

cc <- optimiser(y,x,1000)
optimiser(y,rep(0,t),1000)
optimiser(y,runif(t),1000)
lines(exp(cc[1,]))

plot(gradient2(cc[1,],y[1,],Kx_inv))

lines(exp(cc[3,]))
plot(exp(f[6,]))


plot(x,exp(cc[6,]))
cc2 <- optimiser(y,x,1000)
#plot(x,y[2,])
plot(x,exp(cc[2,]))
plot(x,y[2,])
plot(x,exp(f[2,]))
Kx <- fillup(x,x,kern,1)

X <- x
Kx <- fillup(X,X,kern,1)
L <- t(chol(Kx+diag(t)*10e-8))
#Linv <- backsolve(L, diag(dim(L)[1]))
Linv <- solve(L)
Kx_inv <- t(Linv)%*%Linv
