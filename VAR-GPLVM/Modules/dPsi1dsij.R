dPsi1dsij <- function(i,j,mod){
  t <- dim(mod$mu)[1]
  m <- dim(mod$xu)[1]
  
  end <- matrix(0,ncol=m,nrow=t)
  factor <- numeric(m)
  for (k in 1:m){
    factor[k] <- 0.5*((mod$thetaf$lengthscales[j]*(mod$mu[i,j]-mod$xu[k,j]))^2-0.5*mod$thetaf$lengthscales[j]/(mod$thetaf$lengthscales[j]*mod$Lambda[i,j]+1))
  }
  end[i,] <- factor
  return(mod$psi1*end)
}
