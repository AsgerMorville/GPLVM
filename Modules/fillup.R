#Fillup function
fillup <- function(x1,x2,kernel,theta){
  n1 <- length(x1)
  n2 <- length(x2)
  matr2 <- matrix(NA,nrow=n1,ncol=n2)
  for (i in 1:n1){
    for (j in 1:n2){
      matr2[i,j] <- kernel(x1[i],x2[j],theta)
    }
  }
  return(matr2)
}
