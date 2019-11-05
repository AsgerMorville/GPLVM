#Number of neurons
n <- 10

t <- 100

#Path of mouse
x <- density(x=c(10,90),bw=10,n=t)
x <- x$y
x <- x*100/(0.8)-1.5
plot(x)
phis <- runif(n,max=2*pi)
plot(phis)
omegas <- 1+runif(n,max=3)
f <- matrix(NA,nrow=n,ncol=t)
for (i in 1:n){
  f[i,] <- sin(omegas[i]*x+phis[i])
}
plot((h[8,]))
h <- exp(f)

observed <- f
for (i in 1:n){
  observed[i,] <- rpois(t,lambda=h[i,])
}
plot(observed[1,])
lines(h[1,])

plot(observed[2,])
lines(h[2,])
plot(x,h[5,])

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
  first <- observed[1,]-exp(f)
  second <- invx%*%f
  return(first-second)
}
gradascent <- function(grad,init,gamma){
  a <- init
  for (i in 1:5){
    
    nab <- grad(a)
    a <- a+gamma*nab
    #if (norm(grad(a),type="2")<0.01) {break}
  }
  return(a)
}

gradient(f[1,])
fmax<-gradascent(gradient,startg,0.0001)

A <- exp(diag(fmax))+invx

plot(fmax)
plot(xq)
plot(x,f[2,])
