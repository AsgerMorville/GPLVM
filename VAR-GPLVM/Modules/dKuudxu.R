dKuudxu <- function(i,j,mod){
  #difkernq reimplementation
  #difkernq <- function(index,j,x,K,l){
  Kuu <- mod$Kuu
  l <- mod$thetaf$lengthscales
  x <- mod$xu
  n <- dim(x)[1]
  x <- x[,j]
  output <- matrix(0,nrow=n,ncol=n)
  vec <- x[i]-x
  output[i,] <- Kuu[i,]
  output[,i] <- Kuu[,i]
  output[i,] <- output[i,]*(x-x[i])/l^2
  output[,i] <- output[,i]*(x-x[i])/l^2
  
  output[i,i] <- 0
  return(output)
}