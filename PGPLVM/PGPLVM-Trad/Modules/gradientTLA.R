gradientTLA <- function(model2){
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
  
  
  for (p in 1:N){
    for (j in 1:q){
      Wp <- diag(exp(fhat[,p]))
      inv <- solve(diag(t)+Kx%*%Wp)
      
      
      
      one <- solve(Kxinv+Wp)
      for (i in 1:t){
        dKx <- difkernq(i,j,X,Kx,1)
        first <- 0.5*fhat[,p]%*%Kxinv%*%dKx%*%Kxinv%*%fhat[,p]
        prod <- dKx%*%Wp
        second <- 0.5*sum(diag(inv%*%prod))
        #output[i] <- output[i]+ first-second
        
        #Calculate dlogqdfpi (FOR TLA)
        
        dLogqDf <- -0.5*one[i,i]*(-exp(fhat[i,p]))
        #calkculate dfhatvecdxji
        nabla <- as.vector(y[,p]-exp(fhat[,p]))
        dfhatdx <- inv%*%dKx%*%nabla
        TLA <- dLogqDf*dfhatdx[i]
        output[i,j] <- output[i,j] + first-second+TLA
      }
    }
    
    
  }
  
  return((output)-Ktinv%*%X)
}
