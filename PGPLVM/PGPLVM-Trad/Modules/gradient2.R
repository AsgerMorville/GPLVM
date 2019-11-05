gradient2 <- function(model2){
  #Returns the gradient. Output format: matrix. Input: x = matrix.
  X <- model2$X 
  y <- model2$y
  fhat <- model2$fhat 
  Kx <- model2$Kx 
  Kxinv <- model2$Kx_inv 
  Ktinv <- model2$Ktinv 
  N <- dim(y)[2]
  t <- dim(y)[1]
  q <- dim(X)[2]
  
  output <- matrix(0,nrow=t,ncol=q)
  
  for (i in 1:t){
    for (j in 1:q){
      dKx <- difkernq(i,j,X,Kx,1)
      for (p in 1:N){
        first <- 0.5*fhat[,p]%*%Kxinv%*%dKx%*%Kxinv%*%fhat[,p]
        Wj <- diag(exp(fhat[,p]))
        inv <- solve(diag(t)+Kx%*%Wj)
        prod <- dKx%*%Wj
        second <- 0.5*sum(diag(inv%*%prod))
        output[i] <- output[i]+ first-second
      }
    }
    
    
  }
  
  return(output-Ktinv%*%X)
}
