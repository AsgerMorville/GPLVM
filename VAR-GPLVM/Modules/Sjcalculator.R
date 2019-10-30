
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
  #Calculate S_arr,S_mat
  Kx <- fillup2(tvec,kern=matern,theta=modeltest$thetax)
  UKx <- chol(Kx)
  Ukxinv <- backsolve(r=UKx,x=diag(ncol(UKx)))
  Kxinv <- Ukxinv%*%t(Ukxinv)
  S_arr <- array(NA,dim=c(n,n,q))
  S_mat <- matrix(NA,nrow=n,ncol=q)
  
  for (i in 1:q){
    U1 <- chol(Kxinv+diag(as.vector(L[,i])))
    U1inv <- backsolve(r=U1,x=diag(ncol(U1)))
    S_arr[,,i] <-  U1inv%*%t(U1inv)
    S_mat[,i] <- diag(as.matrix(S_arr[,,i]))
    #Sdets[i] <- 2*sum(log(diag(U1inv))) 
  }
  return(list(S_arr=S_arr,S_mat=S_mat))
}