n <-100

x1 <- runif(n/2)
x2 <- runif(n/2,min=9,max=10)
x <- rnorm(100,mean=seq(from=-5,to=5,length.out = 100),sd=0.5)
plot(x)
x <- c(x1,x2)
plot(x)
kernen <- fillup(x,x)


y1 <- mvrnorm(n = 1, mu=rep(0,n), Sigma=kernen, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
y2 <-mvrnorm(n = 1, mu=rep(0,n), Sigma=kernen, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)


plot(y1,y2)

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
