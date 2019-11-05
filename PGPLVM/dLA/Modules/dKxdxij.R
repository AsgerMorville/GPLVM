dKxdxij <- function(i,j,K,thetaf,xmat){
  tt <- dim(xmat)[1]
  index <- i
  x <- xmat[,j]
  output <- matrix(0,nrow=tt,ncol=tt)
  vec <- x[index]-x
  output[index,] <- K[i,]
  output[,index] <- K[,i]
  output[index,] <- output[index,]*(x-x[index])/thetaf$delta^2
  output[,index] <- output[,index]*(x-x[index])/thetaf$delta^2
  
  return(output)
}


