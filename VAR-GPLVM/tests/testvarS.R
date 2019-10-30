#matrix tester
#This function output is a matrix -> make into vector and apply numderiv

library(numDeriv)
#Function to be tested:
#Test s-derivs
#Our guess:
ourgrad <- function(s){
  mod <- list()
  mod$thetaf$sigmaf <- 1
  S <- as.matrix(rep(1,30))
  S[16,1] <- s
  mod$S_mat <- S
  mod$xu <- data$x
  #mod$xu[3,1] <- xus
  mod$thetaf$lengthscales <-0.1
  mod$mu <- mustart
  
  mat <- dPsi2dsij(16,1,mod)
  
  return(as.vector(mat))
}
ourgrad(2)

testfunc <- function(s){
  thetaf=list(sigmaf=1,lengthscales=0.1)
  S <- as.matrix(rep(1,30))
  S[16,1] <- s
  mu <- mustart
  #modeltest$Lambda <- Lstart
  #modeltest$Lambda[3,1] <- s
  xu <- data$x
  n <- 30
  m <- 30
  #mod1 <- modCalculator(modeltest)
  #mu <- mod1$mu
  #lambda <- Lstart
  return(as.vector(psistats(mu,S,xu,thetaf,n,m)$psi2))
}

numer <- jacobian(testfunc,10)
reel <- ourgrad(10)
plot(numer-reel)
#Endtest <- matrix(numer,nrow=30,ncol=30)

#Endtest

ourgradd <- ourgrad(10)
#endOur <- matrix(ourgradd,nrow=30,ncol=30)
#Endtest-endOur
plot(numer-ourgradd)
numer-ourgradd
#Our dpsi2dxu is okay

#How about dpsi1dxu
#Our guess:
ourgrad <- function(s){
  mod <- list()
  mod$thetaf <- list(sigmaf=1,lengthscales=0.1)
  S <- as.matrix(rep(1,30))
  S[14,1] <- s
  mod$S_mat <- S
  mod$xu <- data$x
  mod$mu <- mustart
  mod$psi1 <- psistats(mod$mu,S,mod$xu,mod$thetaf,30,30)$psi1
  
  mat <- dPsi1dsij(14,1,mod)
  
  return(as.vector(mat))
}

#matrix(ourgrad(2),nrow=30,ncol=30)
testfunc <- function(s){
  thetaf=list(sigmaf=1,lengthscales=0.1)
  S <- as.matrix(rep(1,30))
  S[14,1] <- s
  mu <- mustart
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
numer <- jacobian(testfunc,0.1)
reel <- ourgrad(0.1)
numer
reel
plot(numer-reel)

#dpsi1dxu is also okay

#Lets try dKuudxu
#Our guess
ourgrad <- function(xus){
  mod <- list()
  thetaf=list(sigmaf=1,lengthscales=2)
  
  l <- mod$thetaf$lengthscales[1]
  mod$thetaf <- thetaf
  xu <- data$x
  xu[3,1] <- xus
  mod$xu <- xu
  mod$Kuu <-  fillup(xu,kern=ardkernel,thetaf)
  mat <- dKuudxu(3,1,mod)
  
  return(as.vector(mat))
}
ourgrad(1)

testfunc <- function(xus){
  thetaf=list(sigmaf=1,lengthscales=2)
  xu <- data$x
  xu[3,1] <- xus
  return(as.vector(fillup(xu,kern=ardkernel,thetaf)))
}

#fillup(as.matrix(c(1,2.2,3,4)),kern=ardkernel,thetaf)
reel <- ourgrad(3)
num <- jacobian(testfunc,3)
plot(num-reel)
jacobian(testfunc,1)
plot(ourgrad(10)-jacobian(testfunc,1))
library(numDeriv)
matrix(jacobian(testfunc,1),ncol=30,nrow=30)

#dpsi1dwj
ourgrad <- function(w){
  mod <- list()
  thetaf=list(sigmaf=1,lengthscales=w)
  S <- as.matrix(rep(1,30))
  mod$S_mat <- S
  l <- mod$thetaf$lengthscales[1]
  mod$thetaf <- thetaf
  xu <- data$x
  #xu[3,1] <- xus
  mod$xu <- xu
  mod$mu <- data$x
  mod$y <- data$y
  mod$Kuu <-  fillup(xu,kern=ardkernel,thetaf)
  
  mod$psi1 <- psistats(mod$mu,S,mod$xu,mod$thetaf,30,30)$psi1
  mat <- dPsi1dwj(1,mod)
  return(as.vector(mat))
}

testfunc <- function(w){
  thetaf=list(sigmaf=1,lengthscales=w)
  S <- as.matrix(rep(1,30))
  xu <- data$x
  mu <- data$x
  mat <- psistats(mu,S,xu,thetaf,30,30)$psi1
  return(as.vector(mat))
}

library(numDeriv)
test<- jacobian(testfunc,1)
our <- ourgrad(1)
plot(test-our)

ourgrad <- function(w){
  mod <- list()
  mod$thetaf=list(sigmaf=1,lengthscales=w)
  S <- as.matrix(rep(1,30))
  mod$S_mat <- S
  mod$xu <- data$x
  mod$mu <- data$x
  mod$y <- data$y
  mat <- dPsi2dwj(1,mod)
  
  return(as.vector(mat))
}
ourgrad(1)
testfunc <- function(w){
  thetaf=list(sigmaf=1,lengthscales=w)
  S <- as.matrix(rep(1,30))
  xu <- data$x
  mu <- data$x
  mat <- psistats(mu,S,xu,thetaf,30,30)$psi2
  return(as.vector(mat))
}
num <- jacobian(testfunc,1)
reel <- ourgrad(1)
plot(num-reel)

ourgrad <- function(w){
  mod <- list()
  mod$thetaf=list(sigmaf=1,lengthscales=w)
  mod$xu <-  data$x
  mod$Kuu <-  fillup(mod$xu,kern=ardkernel,mod$thetaf)
  mat <- dKuudwj(1,mod)
  
  return(as.vector(mat))
}
ourgrad(1)
testfunc <- function(w){
  thetaf=list(sigmaf=1,lengthscales=w)
  xu <- data$x
  return(as.vector(fillup(xu,kern=ardkernel,thetaf)))
}
num <- jacobian(testfunc,1)
reel <- ourgrad(1)

plot(num-reel)
###

startcopy <- startguess
#dFhatdS
testfunc <- function(s){
  startguess[40] <-s
  return(vecLogLik(startguess))
}
ourGrad <- function(s){
  #mod <- list()
  
  list <- vecToList(startcopy,t,q,p,m)
  
  modeltest$mubar <- list$mubar
  modeltest$Lambda <- list$Lambda
  modeltest$xu <- list$xu
  modeltest$thetaf$lengthscales <- list$lengthscales
  
  modd <- modCalculator(modeltest)
  modd$
  mod$mu <- as.matrix(startcopy)
  t <- dim(mod$mu)[1]
  q <- dim(mod$mu)[2]
  output <- matrix(NA,nrow=t,ncol=q)
  #Obtain A matrix
  
  Ainv <- mod$Ainv
  YYt <- mod$y%*%t(mod$y)
  sigma <- sqrt(1/(mod$beta))
  psi1 <- mod$psi1
}