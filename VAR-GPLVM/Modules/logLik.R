logLik <- function(model){
  mod <- modCalculator(model)
  
  #Implement sum of (25) and KL term
  Kuu <- mod$Kuu
  sigma <- sqrt(1/mod$beta)
  psi0 <- mod$psi0
  psi1 <- mod$psi1
  psi2 <- mod$psi2
  Kuuinv <- mod$Kuuinv
  thetax <- mod$thetax
  L <- mod$Lambda
  y <- model$y
  tvec <- model$tvec
  mu <- mod$mu
  
  #Calculate W
  L2 <- t(chol((1/sigma^2)*psi2+Kuu+0.001*diag(n)))
  L2inv <- solve(L2)
  L2inv <- t(L2inv)%*%L2inv
  W <- (1/sigma^2)*diag(n)-(1/sigma^4)*psi1%*%L2inv%*%t(psi1)
  
  #quadratic form sum
  sum <- 0
  for (i in 1:p){
    sum <- sum+t(y[,i])%*%W%*%y[,i]
  }
  sum <- -0.5*sum
  #ready to calculate final log likehood expression
  loglik <- -p*n*log(sigma)+(p/2)*2*sum(log(diag(mod$KuuL)))-(p/2)*2*sum(log(diag(L2)))+sum
  
  #add terms outside of log in (25)
  loglik <- loglik-p*psi0/(2*sigma^2)+p*1/(2*sigma^2)*sum(diag(Kuuinv%*%psi2))
  
  #Calculate KL divergence term. equation from page 17
  #First we need Kx, which we call Kt
  Kt <- fillup(tvec,kern=matern,thetax)
  L3 <- t(chol(Kt))
  Ktinv <- solve(L3)
  Ktinv <- t(Ktinv)%*%Ktinv
  
  #Calculate S-matrix from the Lambda matrix
  S <- array(NA,dim=c(n,n,q))
  for (i in 1:q){
    S[,,i] <- solve(Ktinv+diag(as.vector(L[,i])))
  }
  
  KL <- 0
  for (j in 1:q){
    Sj <- as.matrix(S[,,j])
    Sjdet <- as.numeric(determinant(Sj,logarithm=T)$modulus)
    muj <- as.vector(mu[,j])
    KL <- KL +sum(diag(Ktinv%*%Sj+Ktinv%*%mu[,j]%*%t(mu[,j])))+2*sum(log(diag(L3)))-Sjdet
  }
  return(loglik-0.5*KL)
}
