dFdthetaf <- function(mod){
  #Only learn the lengthscales. This vector will be q-long
  #Out put is a q-long vector of derivative wrt w_j for j 1 to q
  q <- dim(mod$mubar)[2]
  output <- numeric(q)
  
  #Precalculate trace terms
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
  
  #Loop over the q-entries
  for (j in 1:q){
    dPsi1dwj <- dPsi1dwj(j,mod)
    dPsi2dwj <- dPsi2dwj(j,mod)
    dKuudwj <- dKuudwj(j,mod)
    psi1term <- mod$beta*sum(diag(t(dPsi1dwj)%*%psi1trace))
    kuuterm <- 0.5*sum(diag(dKuudwj%*%kuutrace))
    psi2term <- 0.5*mod$beta*sum(diag(dPsi2dwj%*%psi2trace))
    output[j] <- psi1term+kuuterm+psi2term
  }
  return(output)
}
