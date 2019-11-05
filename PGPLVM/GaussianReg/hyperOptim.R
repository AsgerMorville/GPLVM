#optimize over kernel hyper parameter
# kernel: 
repeatt <- 100
n <- 50
x <- seq(from=-5,to=5,length.out = n)
#x <- runif(n,min=-5,max=5)
kern <- function(x1,x2,theta){
  exp(-(x1-x2)^2/(2*theta))
}
difkern <-  kern(x,x,theta=1)
difKern2 <- function(x1,x2,theta){
    exp(-(x1-x2)^2/(2*theta))*(x1-x2)^2/(2*theta^2)
}
fillup <- function(x1,x2,kernel,theta){
  n1 <- length(x1)
  n2 <- length(x2)
  matr2 <- matrix(NA,nrow=n1,ncol=n2)
  for (i in 1:n1){
    for (j in 1:n2){
      matr2[i,j] <- kernel(x1[i],x2[j],theta)
    }
  }
  return(matr2)
}




Kernel <- fillup(x,x,kern,3)
difKern <- fillup(x,x,difKern2,3)
Y <- mvrnorm(n=1,mu=rep(0,n),Sigma = Kernel)
Y2 <- Y+rnorm(n,sd=0.1)

plot(x,Y2)
plot(x,Y)
yy <- matrix(nrow=repeatt,ncol=n)
###
for (i in 1:repeatt){
  Y <- mvrnorm(n=1,mu=rep(0,n),Sigma = Kernel)
  Y2 <- Y+rnorm(n,sd=0.1)
  yy[i,] <- Y2
}


gradientHyp <- function(theta,Y22){
  Kernel <- fillup(x,x,kern,theta)
  invx <- solve(Kernel+diag(nrow=n)*0.1^2)
  first <- 0.5*t(Y22)%*%invx
  second <- invx%*%Y22
  Firstterm <- first%*%difKern%*%second
  Second <- 0
  for (i in 1:n){
    Second <- Second + invx[i,]%*%difKern[,i]
  }
  Second <- Second*(-0.5)
  return(Firstterm+Second)
}

for ( i in 1:repeatt){
  
}
source("functions/GAscent.R")
n2 <- 60

x2<- seq(from=1,to=20,length.out=n2)
val <- numeric(n2)
val
for ( i in 1:n2){
    val[i] <-gradientHyp(x2[i],Y2)
  
  print(i)
  
}
val
x2[7]
plot(val)
GAscent(gradientHyp,init=25,0.1)
Checker <- function(theta){
  Kf <- fillup(x,x,kern,theta)
  Ky <- Kf+diag(nrow=n)*0.1^2
  invx <- solve(Ky)
  
}

x2[15]
init <- 1

