fillup <- function(x,kern=ardkernel,theta){
  #x matrix, kern function, thetea 3-long vec
  q<- dim(x)[2]
  n <- dim(x)[1]
  mat <- array(NA,dim=c(n,n,q))
  for (i in 1:q){
    mat[,,i] <- matrix(x[,i],n,n,byrow=F)-matrix(x[,i],n,n,byrow=T)
  }
  mat <- apply(mat,c(1,2),kern,theta)
  #mat <- apply(mat,c(1,2),sum)
  #fin <- apply(mat,c(1,2),kern,theta)
  return(mat)
}

#Fillup2 for matern kernel
fillup2 <- function(x,kern=ardkernel,theta){
  #x matrix, kern function, thetea 3-long vec
  q<- dim(x)[2]
  n <- dim(x)[1]
  mat <- array(NA,dim=c(n,n,q))
  for (i in 1:q){
    mat[,,i] <- abs((matrix(x[,i],n,n,byrow=F)-matrix(x[,i],n,n,byrow=T)))
  }
  mat <- apply(mat,c(1,2),kern,theta)
  #mat <- apply(mat,c(1,2),sum)
  #fin <- apply(mat,c(1,2),kern,theta)
  return(mat)
}
