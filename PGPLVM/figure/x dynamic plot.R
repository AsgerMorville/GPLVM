library(MASS)
library(ggplot2)
#Create example
#F is nonlinear, eg. sin
n <- 20
t <- 100
t2 <- 300
tvec <- seq(from=1,to=100,length.out = 300)
matern <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}
Kt <- fillup(tvec,tvec,matern,20)
#x <- matrix(NA,nrow=1,ncol=t)
x<-mvrnorm(n=1,mu=rep(0,t2),Sigma=Kt)

aa <- qplot(tvec,x,geom = "line")
ggsave(filename="quickplot",aa,device = "pdf")
