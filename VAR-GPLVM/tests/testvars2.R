#Testfunc is 
testfunc <- function(s){
  #mu,Lambda,xu,thetaf,n,m
  mu <- data$x
  Lambda <- matrix(1,ncol=1,nrow=30)
  Lambda[10,1] <- s
  xu <- data$x
  n <- 30
  m <- 30
  return(as.vector(psistats(mu,Lambda,xu,thetaf,n,m)$psi2))
}
testfunc(1)

#Our guess
ourgrad <- function(s){
  #wj <- mod$thetaf$lengthscales[j]

  
  mod<- list()
  mod$xu <- data$x
  mod$mu <- data$x
  mod$thetaf <- thetaf
  mod$S_mat <- matrix(1,ncol=1,nrow=30)
  mod$S_mat[10,1] <- s
  return(dPsi2dsij(10,1,mod))
}

ourgrad(2)
res <- numeric(30)
for (i in 1:30){
  #Testfunc is 
  testfunc <- function(s){
    #mu,Lambda,xu,thetaf,n,m
    mu <- data$x
    Lambda <- matrix(1,ncol=1,nrow=30)
    Lambda[i,1] <- s
    xu <- data$x
    n <- 30
    m <- 30
    return(as.vector(psistats(mu,Lambda,xu,thetaf,n,m)$psi2))
  }
  testfunc(1)
  
  #Our guess
  ourgrad <- function(s){
    #wj <- mod$thetaf$lengthscales[j]
    
    
    mod<- list()
    mod$xu <- data$x
    mod$mu <- data$x
    mod$thetaf <- thetaf
    mod$S_mat <- matrix(1,ncol=1,nrow=30)
    mod$S_mat[i,1] <- s
    return(dPsi2dsij(i,1,mod))
  }
  num <- jacobian(testfunc,1)
  reel <- as.vector(ourgrad(1))
  
  res[i] <- max(abs(num-reel))
  print(i)
}
num <- jacobian(testfunc,0.01)
reel <- as.vector(ourgrad(0.01))
plot(num-reel)
plot(res)
##How about dPsi1
testfunc <- function(s){
  #mu,Lambda,xu,thetaf,n,m
  mu <- data$x
  Lambda <- matrix(1,ncol=1,nrow=30)
  Lambda[14,1] <- s
  xu <- data$x
  n <- 30
  m <- 30
  return(as.vector(psistats(mu,Lambda,xu,thetaf,n,m)$psi1))
}
testfunc2 <- function(s){
  thetaff=list(sigmaf=1,lengthscales=1)
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
  return(as.vector(psistats(mu,S,xu,thetaff,n,m)$psi1))
}

ourgrad <- function(s){
  #wj <- mod$thetaf$lengthscales[j]
  
  mod<- list()
  mod$xu <- data$x
  mod$mu <- data$x
  mod$thetaf <- thetaf
  mod$S_mat <- matrix(1,ncol=1,nrow=30)
  mod$S_mat[14,1] <- s
  mod$psi1 <- psistats(mod$mu,mod$S_mat,mod$xu,thetaf,30,30)$psi1
  return(dPsi1dsij(14,1,mod))
}
testfunc2(1)
num<-jacobian(testfunc,0.1)
rel<- ourgrad(0.1)
num
rel
plot(num-as.vector(rel))
