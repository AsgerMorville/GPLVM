eigenInit <- function(y,latentDim){
  YYt <- y%*%t(y)
  q <- latentDim
  alpha <- 1
  beta <- 1
  n <- dim(y)[1]
  eig <- eigen(YYt)
  U <- matrix(NA,nrow=n,ncol=q)
  for (i in 1:q){
    U[,i] <- eig$vectors[,i]
  }
  
  D <- dim(y)[2]
  l <- matrix(0,ncol=q,nrow=q)
  for (i in 1:q){
    l[i,i] <- 1/sqrt(eig$values[i]/alpha*D-1/alpha*beta)
  }
  return(U%*%l)
}
