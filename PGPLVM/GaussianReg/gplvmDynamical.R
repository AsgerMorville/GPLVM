library(MASS)
library(ggplot2)


plot(density(x=c(1,2),from=0,to=10))
x1 <-density(x=c(1,2),from=0,to=10,n=100)

x1 <- x1$y

x2 <-density(x=c(8,9),from=0,to=10,n = 100)
plot(x2)
x2 <- x2$y
plot(1:100, x1, type="b",col="red" )
par(new=TRUE)
plot(1:100, x2, col="green" ,type="b")

kern <- function(v1,v2){
  vec <- v1-v2
  return(0.1*exp(-(norm(vec,"2"))^2/20))
}

leng <- 100
matr2 <- matrix(NA,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    matr2[i,j] <- kern(c(x1[i],x2[i]),c(x1[j],x2[j]))
  }
}
h <-mvrnorm(n = 1, mu=rep(0,100), Sigma=matr2)
par(new=TRUE)
plot(h,type="b")
