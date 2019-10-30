#Calculate derivative matrix of Fhat wrt mumean

dFhatdmu <- function(mod){
  #Output is matrix of size t times q
  
  t <- dim(mod$mu)[1]
  q <- dim(mod$mu)[2]
  output <- matrix(NA,nrow=t,ncol=q)
  #Obtain A matrix
  sigma <- sqrt(1/(mod$beta))
  psi1 <- mod$psi1
  Ainv <- mod$Ainv
  YYt <- mod$y%*%t(mod$y)
  
  #calculate expressions inside trace terms as they dont change with i,j
  psi1trace <- YYt%*%psi1%*%Ainv
  psi2trace <- p*mod$Kuuinv-sigma^2*p*Ainv-Ainv%*%t(psi1)%*%YYt%*%psi1%*%Ainv
  
  for (i in 1:t){
    for (j in 1:q){
      #Note that first term is zero
      dPsi1dmuj <- dPsi1dmu(i,j,mod)
      dPsi2dmuj <- dPsi2dmu(i,j,mod)
      second <- mod$beta*sum(diag(t(dPsi1dmuj)%*%psi1trace))
      third <- 0.5*mod$beta*sum(diag(dPsi2dmuj%*%psi2trace))
      output[i,j] <- second+third
    }
  }
  return(output)
}
