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
  L2 <- t(chol((1/sigma^2)*psi2+Kuu))
  U22 <- t(L2)
  U22inv <- backsolve(r=U22,x=diag(ncol(U22)))
  Wintinv <- U22inv%*%t(U22inv)
  
  W <- (1/sigma^2)*diag(n)-(1/sigma^4)*psi1%*%Wintinv%*%t(psi1)
  
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
  Kt <- fillup2(tvec,kern=matern,thetax)
  L3 <- t(chol(Kt))
  Ktinv <- solve(L3)
  Ktinv <- t(Ktinv)%*%Ktinv
  
  #Calculate S-matrix from the Lambda matrix
  S <- array(NA,dim=c(n,n,q))
  Sdets <- numeric(q)
  for (i in 1:q){
    U1 <- chol(Ktinv+diag(as.vector(L[,i])))
    U1inv <- backsolve(r=U1,x=diag(ncol(U1)))
    S[,,i] <-  U1inv%*%t(U1inv)
    #print(diag(U1inv))
    #print(2*sum(log(diag(U1inv))))
    Sdets[i] <- 2*sum(log(diag(U1inv))) 
  }
  Sdets1 <- numeric(q)
  KL <- 0
  for (j in 1:q){
    Sj <- as.matrix(S[,,j])
    Sjdet <-Sdets[j]
    muj <- as.vector(mu[,j])
    KL <- KL +sum(diag(Ktinv%*%Sj+Ktinv%*%mu[,j]%*%t(mu[,j])))+2*sum(log(diag(L3)))-Sjdet
  }
  
  return(loglik-0.5*KL)
}

