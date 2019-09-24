modCalculator <- function(model){
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
  
  sigma <- 1/sqrt(abs(beta))
  sigmaf <- thetaf$sigmaf
  

  #calculate Kx
  Kx <- fillup(tvec,kern=matern,theta=model$thetax)
  
  #calculate variational means
  mu <- Kx%*%mubar
  
  #calculate psi-statistics
  Psistats <- psistats(mu,Lambda,xu,thetaf,n,m)
  psi0 <- Psistats$psi0
  psi1 <- Psistats$psi1
  psi2 <- Psistats$psi2
  
  #Find Kuu, and Lower cholesky x2
  #theta <- list(sigmaf=sigmaf,lengthscales=w)
  Kuu <- fillup(xu,kern=ardkernel,thetaf)+0.0001*diag(m)
  L <- t(chol(Kuu))
  #L2 <- t(chol((1/sigma^2)*psi2+Kuu))
  
  Kuuinv <- solve(L)
  Kuuinv <- t(Kuuinv)%*%Kuuinv
  
  mod <- list()
  #Assign the variables to mod
  mod$tvec <- tvec
  mod$y <- y
  mod$mubar <- mubar
  mod$mu <- mu
  mod$Lambda<- Lambda 
  mod$xu <- xu
  mod$psi0 <- psi0
  mod$psi1 <- psi1
  mod$psi2 <- psi2
  mod$KuuL <- L
  mod$Kuu <- Kuu
  mod$Kuuinv <- Kuuinv
  mod$beta <- beta
  mod$thetaf <- thetaf
  mod$thetax <- thetax
  mod$Kx <- Kx
  
  
  return(mod)
}
