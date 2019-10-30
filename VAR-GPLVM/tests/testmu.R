ourgrad <- function(mus){
  mod <- list()
  mod$thetaf$sigmaf <- 1
  S <- as.matrix(rep(1,30))
  mod$S_mat <- S
  mod$xu <- data$x
  #mod$xu[3,1] <- xus
  mod$thetaf$lengthscales <-2
  mod$mu <- mustart
  mod$mu[3,1] <- mus
  
  mat <- dPsi2dmu(3,1,mod)
  
  return(as.vector(mat))
}

testfunc <- function(mus){
  thetaf=list(sigmaf=1,lengthscales=2)
  S <- as.matrix(rep(1,30))
  mu <- mustart
  mu[3,1] <- mus
  #modeltest$Lambda <- Lstart
  #modeltest$Lambda[3,1] <- s
  xu <- data$x
  #xu[3,1] <- xus
  n <- 30
  m <- 30
  #mod1 <- modCalculator(modeltest)
  #mu <- mod1$mu
  #lambda <- Lstart
  return(as.vector(psistats(mu,S,xu,thetaf,n,m)$psi2))
}

numer <- jacobian(testfunc,10)
#Endtest <- matrix(numer,nrow=30,ncol=30)

#Endtest

ourgradd <- ourgrad(10)
#endOur <- matrix(ourgradd,nrow=30,ncol=30)
#Endtest-endOur
plot(numer-ourgradd)
#
ourgrad <- function(mus){
  mod <- list()
  mod$thetaf$sigmaf <- 1
  S <- as.matrix(rep(1,30))
  mod$S_mat <- S
  mod$xu <- data$x
  #mod$xu[3,1] <- xus
  mod$thetaf$lengthscales <-2
  mod$mu <- mustart
  mod$mu[3,1] <- mus
  mod$psi1 <- psistats(mod$mu,S,mod$xu,mod$thetaf,30,30)$psi1
  
  mat <- dPsi1dmu(3,1,mod)
  
  return(as.vector(mat))
}

testfunc <- function(mus){
  thetaf=list(sigmaf=1,lengthscales=2)
  S <- as.matrix(rep(1,30))
  mu <- mustart
  mu[3,1] <- mus
  #modeltest$Lambda <- Lstart
  #modeltest$Lambda[3,1] <- s
  xu <- data$x
  #xu[3,1] <- xus
  n <- 30
  m <- 30
  #mod1 <- modCalculator(modeltest)
  #mu <- mod1$mu
  #lambda <- Lstart
  return(as.vector(psistats(mu,S,xu,thetaf,n,m)$psi1))
}

testfunc(1)
numer <- jacobian(testfunc,1)
Endtest <- matrix(numer,nrow=30,ncol=30)
Endtest-matrix(ourgrad(1),nrow=30,ncol=30)

plot(numer-ourgrad(1))
