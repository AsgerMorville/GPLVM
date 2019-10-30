
Sjcalculator <- function(L,tvec,n,q){
  Kx <- fillup2(tvec,kern=matern,theta=modeltest$thetax)
  Kxinv <- solve(Kx)
  output <- array(0,dim=c(n,n,q))
  output2 <- matrix(NA,nrow=n,ncol=q)
  for (j in 1:q){
    output[,,j] <- solve(as.matrix(Kxinv+diag(as.vector(L[,j]))))
    output2[,j] <- pmax(diag(as.matrix(output[,,j])),0.05)
  }
  #Extract diagonals
  
  return(as.matrix(output2))
}

Sjcalculator2 <- function(L,tvec,n,q){
  Kx <- fillup2(tvec,kern=matern,theta=modeltest$thetax)
  Kxinv <- solve(Kx)
  output <- array(0,dim=c(n,n,q))
  output2 <- matrix(NA,nrow=n,ncol=q)
  for (j in 1:q){
    output[,,j] <- solve(as.matrix(Kxinv+diag(as.vector(L[,j]))))
    output2[,j] <- diag(as.matrix(output[,,j]))
  }
  #Extract diagonals
  
  return(as.matrix(output2))
}