#dPsi1dwj

dPsi1dwj <- function(j,mod){
  psi1 <- mod$psi1
  m <- dim(mod$xu)[1]
  t <- dim(mod$y)[1]
  factor <- matrix(NA,ncol=m,nrow=t)
  mu <- mod$mu
  S_mat <- mod$S_mat
  xu <- mod$xu
  w <- mod$thetaf$lengthscales[j]
  
  for (i in 1:t){
    for (k in 1:m){
      factor[i,k] <- ((mu[i,j]-xu[k,j])/(w*S_mat[i,j]+1))^2+S_mat[i,j]/(w*S_mat[i,j]+1)
    }
  }
  factor <- factor*(-0.5)
  return(psi1*factor)
  
}
