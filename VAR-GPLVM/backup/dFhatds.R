dFhatds <- function(mod){
  #output: returns the matrix of derivative dFhat/ds
  t <- dim(mod$mu)[1]
  q <- dim(mod$mu)[2]
  output <- matrix(NA,nrow=t,ncol=q)
  #Obtain A matrix
  
  A <- (1/mod$beta)*mod$Kuu+mod$psi2
  Ainv <- solve(A)
  YYt <- mod$y%*%t(mod$y)
  sigma <- sqrt(1/(mod$beta))
  psi1 <- mod$psi1
  #calculate expressions inside trace terms as they dont change with i,j
  psi1trace <- YYt%*%psi1%*%Ainv
  psi2trace <- p*mod$Kuuinv-sigma^2*p*Ainv-Ainv%*%t(psi1)%*%YYt%*%psi1%*%Ainv
  
  for (i in 1:t){
    for (j in 1:q){
      #Note that first term is zero
      dPsi1dsij <- dPsi1dsij(i,j,mod)
      dPsi2dsij <- dPsi2dsij(i,j,mod)
      second <- mod$beta*sum(diag(t(dPsi1dsij)%*%psi1trace))
      third <- 0.5*mod$beta*sum(diag(dPsi2dsij%*%psi2trace))
      output[i,j] <- second+third
    }
  }
  return(output)
}
