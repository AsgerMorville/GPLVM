#Gaus process
n <- 100
xgrid <- runif(n,min=-5,max=5)

kern <- function(x1,x2){
  exp(-(x1-x2)^2/2)
}

matr <- matrix(NA,nrow=n,ncol=n)

for (i in 1:n){
  for (j in 1:n){
    matr[i,j] <- kern(xgrid[i],xgrid[j])
  }
}
ys <- mvrnorm(n = 1, mu=rep(0,n), Sigma=matr, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#add noise
y1 <- rnorm(n,mean=ys,sd=0.3)

plot(xgrid,ys)


#Gaussian regression
newx <- runif(n,min=-5,max=5)


matr2 <- fillup(xgrid,xgrid)

s2 <- 0.3^2
eye <- diag(x=s2,nrow=n,ncol=n)
mellem <- solve(matr2+eye)

est <- matr2%*%mellem
est2 <- est%*%y1




plot( x, y2, type="l", col="green" )
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
kss <- fillup(newx,newx)
ksx <- fillup(newx,xgrid)
kxs <- t(ksx)
kxx <- fillup(xgrid,xgrid)
covfin <- kss-(ksx%*%mellem%*%kxs)

samp <-  mvrnorm(n = 1, mu=est2, Sigma=covfin, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
plot(xgrid,samp)
## NOW USE ONLY 10 POINTS, evenly spaced
sigm <- 0.3
n2 <- 40
newgrid <- runif(n2,min=-5,max=5)
Kmm <- fillup(newgrid,newgrid)
Kmn <- fillup(newgrid,xgrid)
Knm <- t(Kmn)
summ <- Kmm+sigm^(-2)*Kmn%*%Knm
sigmMat <- solve(summ+diag(x=1e-08,nrow = n2,ncol=n2))
muu <- sigm^(-2)*Kmm%*%sigmMat%*%Kmn%*%y1
plot(newgrid,muu)
par(new=TRUE)
plot(xgrid,est2)
