#MAp estimate of f
x <- 1:100
f <-(density(x=5,from = 0,to=10,bw=1,n=100))
f <- f$y*5
plot(exp(f))
pointss <- matrix(NA,nrow=100,ncol=100)


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
plot(x)
kx <- fillup(x,x)
invx <-solve(kx)

for (j in 1:100){
  observed <- rpois(100,lambda=exp(f))
  y <- observed
  gradient <- function(f){
    first <- y-exp(f)
    second <- invx%*%f
    return(first-second)
  }
  gradascent <- function(grad,init,gamma){
    a <- init
    for (i in 1:500){
      nab <- grad(a)
      a <- a+gamma*nab
      if (norm(grad(a),type="2")<0.1) {break}
    }
    return(a)
  }
  pointss[j,]<-gradascent(gradient,runif(100,min=0.2,max=5),0.01)
  print(j)
}



end <- apply(pointss,c(2),mean)
plot(end)
lines((f))
#


x <- 1:100



plot(f)
plot(observed)
lines(f)


plot(fin)
norm(gradient(fin),type="2")
