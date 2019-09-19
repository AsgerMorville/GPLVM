fillup2 <- function(x1,kern,theta){
  n <- dim(x1)[1]
  q <- dim(x1)[2]
  mat <- array(NA,dim=c(n,n,q))
  for (i in 1:q){
    mat[,,i] <- (matrix(x1[,i],n,n,byrow=F)-matrix(x1[,i],n,n,byrow=T))^2
  }
  mat <- apply(mat,c(1,2),sum)
  fin <- apply(mat,c(1,2),kern,theta)
  return(fin)
}
