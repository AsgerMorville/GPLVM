dFdlambda <- function(mod){
  #calculates the matrix of derivatives of Lambda matrix-function
  Kxinv <- solve(mod$Kx)
  n <- dim(mod$mu)[1]
  q <- dim(mod$mu)[2]
  output <- matrix(NA,ncol=q,nrow=n)
  
  dFhatds <- dFhatds(mod)
  for (j in 1:q){
    #Calculate Sj
    Sj <- solve(Kxinv+diag(as.vector(mod$Lambda[,j])))
    
    #Calculate dFhatdsj
    dFhatdsj <- dFhatds[,j]
    
    #Obtain the j'th vector
    output[,j] <- -(Sj%*%Sj)%*%(dFhatdsj+0.5*mod$Lambda[,j])
  }
  return(output)
}