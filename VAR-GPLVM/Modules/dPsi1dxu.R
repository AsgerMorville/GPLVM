dPsi1dxu <- function(i,j,mod){
  m <- dim(mod$xu)[1]
  t <- dim(mod$mu)[1]
  output <- matrix(0,nrow=t,ncol=m)
  
  #Find the muvector,w,S vector
  w <- mod$thetaf$lengthscales[j]
  mui <- as.vector(mod$mu[,j])
  Si <- as.vector(mod$Lambda[,j])
  psi1 <- mod$psi1
  xumq <- as.numeric(mod$xu[i,j])
  factor_1 <- w*(mui-xumq)/(w*Si+1)
  output[,j] <- factor_1
  return(psi1*output)
}