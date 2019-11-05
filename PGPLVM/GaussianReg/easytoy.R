#Toy easy X

#Number of neurons
n <- 10

t <- 100

#Path of mouse
x <- density(x=c(50),bw=10,n=t)
x <- x$y
x <- x*100/(2)-1.5
plot(x)
phis <- runif(n,max=2*pi)
plot(phis)
omegas <- 1+runif(n,max=3)
f <- matrix(NA,nrow=n,ncol=t)
for (i in 1:n){
  f[i,] <- sin(omegas[i]*x+phis[i])
}
h <- exp(f)
plot((h[6,]))


observed <- f
for (i in 1:n){
  observed[i,] <- rpois(t,lambda=h[i,])
}
plot(observed[1,])
lines(h[1,])

plot(observed[3,])
lines(h[3,])
plot(x,f[3,])

##
startg <- rnorm(t,mean=f[1,],sd=0.1)
plot(startg)
plot(x)
#Step 1. Computer posterior mode and precision matrix
#Gradient ascent
xq <- rnorm(t,mean=x,sd=0)
plot(xq)
invx <- solve(fillup(xq,xq)+diag(t)*10e-6)

fillup <- function(x1,x2){
  n1 <- length(x1)
  n2 <- length(x2)
  matr2 <- matrix(NA,nrow=n1,ncol=n2)
  for (i in 1:n1){
    for (j in 1:n2){
      matr2[i,j] <- kern(x1[i],x2[j])
    }
  }
  return(matr2)
}
kern <- function(x1,x2){
  exp(-(x1-x2)^2/2)
}
gradient <- function(f){
  first <- observed[2,]-exp(f)
  second <- invx%*%f
  return(first-second)
}
gradascent <- function(grad,init,gamma){
  a <- init
  for (i in 1:50){
    
    nab <- grad(a)
    a <- a+gamma*nab
    #if (norm(grad(a),type="2")<0.01) {break}
  }
  return(a)
}

##Newtons method
newton <- function(g,H,start){
  f <- start
  for (i in 1:10000){
    f <- f-0.01*solve(H(f))%*%g(f)
    if (norm(g(f),type="2")<0.1){break}
  }
  return(f)
}
Hessian <- function(f){
  return(diag(exp(f))-invx)
}

fmax<-newton(gradient,Hessian,c(rep(0,100)))
plot(fmax)
plot(f[2,])

checker(f[2,])
##


gradient(f[1,])
fmax<-gradascent(gradient,rep(0,100),0.0001)
plot(observed[1,])
lines(exp(fmax))
A <- exp(diag(fmax))+invx

plot(fmax)
plot(xq)
plot(x,f[2,])

checker <- function(ff){
  return(sum(log(dpois(observed[2,],lambda=exp(ff))))-0.5*t(ff)%*%invx%*%ff)
}
checker(fmax)
checker(rep(0,100))
exp(f[1,])
checker(exp(f[3,]))
plot(observed[1,])
lines(exp(f[3,]))
dpois(observed[1,],lambda=exp(f[1,]))
f[1,]
dpois(observed[1,],lambda=rep(1,100))
checker(observed[1,])
      