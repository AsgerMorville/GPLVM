#Fillup. Has to handle matrices
#fillup <- function(x1,x2,kern,l){
#  n1 <- length(x1)
#  n2 <- length(x2)
#  matr2 <- matrix(NA,nrow=n1,ncol=n2)
#  for (i in 1:n1){
#    for (j in 1:n2){
#      matr2[i,j] <- kern(x1[i],x2[j],l)#
#    }
#  }
#  return(matr2)
#}

fillup <- function(x1,kern,theta){
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