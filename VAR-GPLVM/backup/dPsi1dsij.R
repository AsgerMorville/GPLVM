dPsi1dsij <- function(i,j,mod){
  t <- dim(mod$mu)[1]
  m <- dim(mod$xu)[1]
  w <- mod$thetaf$lengthscales[j]
  muij <- mod$mu[i,j]
  #xuij <- mod$xu[i,j]
  end <- matrix(0,ncol=m,nrow=t)
  Sij <- mod$S_mat[i,j]
  factor <- numeric(m)
  for (k in 1:m){
    factor[k] <- 0.5*(w*(muij-mod$xu[k,j])/(w*Sij+1))^2-0.5*w/(w*Sij+1)
    #factor[k] <- (0.5*(mod$thetaf$lengthscales[j]*(mod$mu[i,j]-mod$xu[k,j]))^2-0.5*mod$thetaf$lengthscales[j]/(mod$thetaf$lengthscales[j]*mod$S_mat[i,j]+1))
  }
  end[i,] <- factor
  return(mod$psi1*end)
}
