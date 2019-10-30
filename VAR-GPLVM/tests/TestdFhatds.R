testFhat <- function(s){
  modeltest$xu <- as.matrix(data$x[indices])
  modeltest$mubar <- data$x
  modeltest$thetaf=list(sigmaf=1,lengthscales=1)
  mod <- modCalc2(modeltest,s)
  
  Kuu <- mod$Kuu
  sigma <- sqrt(1/mod$beta)
  psi0 <- mod$psi0
  psi1 <- mod$psi1
  psi2 <- mod$psi2
  Kuuinv <- mod$Kuuinv
  thetax <- mod$thetax
  L <- mod$Lambda
  y <- data$y
  tvec <- mod$tvec
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
    #U1 <- chol(Ktinv+diag(as.vector(L[,i])))
    U1inv <- chol(mod$S_arr[,,i])
    #S[,,i] <-  U1inv%*%t(U1inv)
    #print(diag(U1inv))
    #print(2*sum(log(diag(U1inv))))
    Sdets[i] <- 2*sum(log(diag(U1inv))) 
  }
  Sdets1 <- numeric(q)
  KL <- 0
  for (j in 1:q){
    Sj <- as.matrix(mod$S_arr[,,j])
    Sjdet <-Sdets[j]
    muj <- as.vector(mu[,j])
    KL <- KL +sum(diag(Ktinv%*%Sj+Ktinv%*%mu[,j]%*%t(mu[,j])))+2*sum(log(diag(L3)))-Sjdet
  }
  
  return(loglik)
}

modCalc2 <- function(model,s){
  tvec <- model$tvec
  y <- model$y
  mubar <- model$mubar
  Lambda <- model$Lambda
  xu <- model$xu
  thetaf <- model$thetaf
  thetax <- model$thetax
  beta <- model$beta
  
  n <- dim(mubar)[1]
  q <- dim(mubar)[2]
  p <- dim(y)[2]
  m <- dim(xu)[1]
  
  sigma <- 1/sqrt(beta)
  sigmaf <- thetaf$sigmaf
  
  
  #calculate Kx
  Kx <- fillup2(tvec,kern=matern,theta=model$thetax)
  
  #calculate variational means
  mu <- Kx%*%mubar
  
  #Calculate S-matrix
  
  S_mat <- objj$S_mat
  S_arr <- objj$S_arr
  S_arr[2,2,1] <- s
  S_mat[2,1] <- s
  #calculate psi-statistics
  Psistats <- psistats(mu,S_mat,xu,thetaf,n,m)
  psi0 <- Psistats$psi0
  psi1 <- Psistats$psi1
  psi2 <- Psistats$psi2
  
  #Find Kuu, and Lower cholesky x2
  #theta <- list(sigmaf=sigmaf,lengthscales=w)
  Kuu <- fillup(xu,kern=ardkernel,thetaf)+diag(m)*100
  U <- chol(Kuu)
  #L2 <- t(chol((1/sigma^2)*psi2+Kuu))
  Uinv <- backsolve(r=U,x=diag(ncol(U)))
  Kuuinv <- Uinv%*%t(Uinv)
  
  #Inv A matrix
  A <- (1/beta)*Kuu+psi2
  #Invert A in a numerically stable way. We try with Choleksy
  LA <- chol(A)
  InvLA <- backsolve(r=LA,x=diag(ncol(LA)))
  Ainv <- (InvLA)%*%t(InvLA)
  
  #Assign the variables to mod
  mod <- list()
  mod$tvec <- tvec
  mod$y <- y
  mod$mubar <- mubar
  mod$mu <- mu
  mod$Lambda<- Lambda 
  mod$xu <- xu
  mod$psi0 <- psi0
  mod$psi1 <- psi1
  mod$psi2 <- psi2
  mod$KuuL <- t(U)
  mod$Kuu <- Kuu
  mod$Kuuinv <- Kuuinv
  mod$beta <- beta
  mod$thetaf <- thetaf
  mod$thetax <- thetax
  mod$Kx <- Kx
  mod$S_mat <- S_mat
  mod$S_arr <- S_arr
  mod$Ainv <- Ainv
  mod$KUUFACTORUP <- Uinv
  return(mod)
}
ourgrad <- function(s){
  modeltest$xu <- as.matrix(data$x[indices])
  modeltest$thetaf=list(sigmaf=1,lengthscales=1)
  #modeltest$thetax=list(l=1)
  modeltest$mubar <- data$x
  mod <- modCalc2(modeltest,s)
  return(dFhatds(mod))
}

num <- grad(testFhat,10)
reel <- ourgrad(10)
num
reel
ourgrad()
