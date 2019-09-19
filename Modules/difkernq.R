#difkernq
difkernq <- function(index,j,x,K,l){
  n <- dim(x)[1]
  x <- x[,j]
  output <- matrix(0,nrow=n,ncol=n)
  vec <- x[index]-x
  output[index,] <- K[index,]
  output[,index] <- K[,index]
  output[index,] <- output[index,]*(x-x[index])/l^2
  output[,index] <- output[,index]*(x-x[index])/l^2
  
  output[index,index] <- 0
  return(output)
}
