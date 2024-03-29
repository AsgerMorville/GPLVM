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
  
  sigma <- 1/sqrt(beta)
  sigmaf <- thetaf$sigmaf
  

  #calculate Kx
  Kx <- fillup2(tvec,kern=matern,theta=model$thetax)
  
  #calculate variational means
  mu <- Kx%*%mubar
  
  #Calculate S-matrix
  S_obj <- Sjcalculator2(Lambda,tvec,n,q)
  S_mat <- S_obj$S_mat
  S_arr <- S_obj$S_arr
  
  #calculate psi-statistics
  Psistats <- psistats(mu,S_mat,xu,thetaf,n,m)
  psi0 <- Psistats$psi0
  psi1 <- Psistats$psi1
  psi2 <- Psistats$psi2
  
  #Find Kuu, and Lower cholesky x2
  #theta <- list(sigmaf=sigmaf,lengthscales=w)
  #Kuu <- fillup(xu,kern=ardkernel,thetaf)+diag(m)*1e-05
  Kuu <- fillup(xu,kern=ardkernel,thetaf)
  #Make loop of tries to factorize kuu matrix.
  for (i in 1:50){
    tryCatch(
      expr = {
        j <- i-4
        U<-chol(Kuu)+diag(m)*10^j
        Kuu <- t(U)%*%U
      },
      error = function(e){ 
        print("adding jiter")
        print(
      }
    )
  }
  
  
  ##
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
