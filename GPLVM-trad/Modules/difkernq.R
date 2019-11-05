#difkernq
difkernq <- function(index,j,x,K,g){
  n <- dim(x)[1]
  x <- x[,j]
  output <- matrix(0,nrow=n,ncol=n)
  vec <- x[index]-x
  output[index,] <- K[index,]
  output[,index] <- K[,index]
  output[index,] <- output[index,]*(x-x[index])*g
  output[,index] <- output[,index]*(x-x[index])*g
  
  output[index,index] <- 0
  return(output)
}
