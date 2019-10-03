dFdxu <- function(mod){
  #returns the matrix (m by q) of derivatives.
  #each poitn consists of sum of three points
  m <- dim(mod$xu)[1]
  q <- dim(mod$mu)[2]
  output <- matrix(NA,nrow=m,ncol=q)
  #Obtain A matrix
  sigma <- sqrt(1/(mod$beta))
  psi1 <- mod$psi1
  A <- (1/mod$beta)*mod$Kuu+mod$psi2
  Ainv <- solve(A)
  YYt <- mod$y%*%t(mod$y)
  Kuuinv <- mod$Kuuinv
  
  #calculate expressions inside trace terms as they dont change with i,j
  psi1trace <- YYt%*%psi1%*%Ainv
  psi2trace <- p*mod$Kuuinv-sigma^2*p*Ainv-Ainv%*%t(psi1)%*%YYt%*%psi1%*%Ainv
  kuutrace <- p*Kuuinv-(1/mod$beta)*p*Ainv-Ainv%*%t(psi1)%*%YYt%*%psi1%*%Ainv-mod$beta*p*Kuuinv%*%mod$psi2%*%Kuuinv
  for (k in 1:m){
    for (j in 1:q){
      #Note that first term is zero
      dPsi1dxu <- dPsi1dxu(k,j,mod)
      dPsi2dxu <- dPsi2dxu(k,j,mod)
      dKuudxu <- dKuudxu(k,j,mod)
      psi1term <- mod$beta*sum(diag(t(dPsi1dxu)%*%psi1trace))
      kuuterm <- 0.5*sum(diag(dKuudxu%*%kuutrace))
      psi2term <- 0.5*mod$beta*sum(diag(dPsi2dxu%*%psi2trace))
      output[k,j] <- psi1term+kuuterm+psi2term
    }
  }
  return(output)
}