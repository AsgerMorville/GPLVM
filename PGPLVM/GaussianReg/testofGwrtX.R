#Create Y's and f's

N <- 20
t <- 100
X <- c(rnorm(t/2,0,sd=0.1),rnorm(t/2,2,sd=0.1))
plot(X)
Ky <- fillup(X,X,kern,1)+diag(t)*0.1^2

#Xn <- density(x=5,from=0,to=10,bw=1,n=t)
#Xn <- Xn$y
#plot(Xn)
yy <- matrix(NA,nrow=N,ncol=t)
for (i in 1:20){
  yy[i,] <- mvrnorm(n=1,mu=rep(0,t),Sigma = Ky)
}

plot(yy[1,])


plot(yy[1,],yy[10,])

#Find initialization
newy <- t(yy)
scnewy <- scale(newy,scale = F)
YYt <- scnewy%*%t(scnewy)
alpha <- 1
beta <- 1

eig <- eigen(YYt)
D <- N
lj <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
Xny <- eig$vectors[,1]*lj

plot(Xny)
###

l <- 1
sigma <- 0.1


loglikchol(c(rep(0,50),rep(9,50)))
loglikchol(c(rep(0,25),c(rep(4,25)),rep(1,50)))
gradientchol(runif(100))

optim(runif(100),loglikchol,gradientchol,method="BFGS")

fine <- 1000
grid <- seq(from=-10,to=20,length.out=fine)
val <- numeric(fine)
for(i in 1:fine){
  val[i]<- loglikchol(c(rep(0,25),rep(0,25),rep(grid[i],50)))
}
plot(grid,val)
grid[66]
loglikchol(runif(100))
library(devtools)
install_github("vqv/ggbiplot")
test <- prcomp(dff)

autoplot(test)
?prcomp


newy <- newy-colMeans(newy)
scnewy <- scale(newy,scale = F)
svd <- svd(scnewy)
col <- svd$u%*%diag(svd$d)

plot(col[,2],col[,1])

d[3:20] <- 0
