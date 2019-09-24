dFdmubar <- function(mod){
  #output a matrix, t times q
  n <- dim(mod$mu)[1]
  q <- dim(mod$mu)[2]
  output <- matrix(NA,nrow=n,ncol=q)
  Kx <- mod$Kx
  
  #Calculate dFhat derivative-matrix
  dFhat <- dFhatdmu(mod)
 
  for (j in 1:q){
    dmuj <- dFhat[,j]
    output[,j] <- Kx%*%(dmuj-mod$mubar[,j])
  }
  return(output)
}
