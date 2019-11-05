#fullLogLik
fullLogLik <- function(model){
  tvec <- model$tvec
  y <- model$y
  X <- model$X
  maxiter <- model$maxiter
  
  N <- dim(y)[2]
  t <- dim(y)[1]
  q <- dim(X)[2]
  
  Ktinv <- model$Ktinv
  #1. a) Computer posterior mode of f_i for i in 1:N
  
  #Compute inverse
  Kx <- fillup(X,ARDKernel,1)
  print(dim(Kx))
  print(t)
  L <- t(chol(Kx+diag(t)*10e-8))
  #Linv <- backsolve(L, diag(dim(L)[1]))
  Linv <- solve(L)
  Kx_inv <- t(Linv)%*%Linv
  
  fhat <- matrix(NA,nrow=t,ncol=N)
  for (i in 1:N){
    fhat[,i] <- newton(y[,i],X,gradientNewton,Hessian,maxiter,Kx_inv)
    #print(i)
  }
  print("done")
  sum <- 0
  for (i in 1:N){
    one <- sum(y[,i]*fhat[,i])-sum(exp(fhat[,i]))
    two <- 0.5*t(fhat[,i])%*%Kx_inv%*%(fhat[,i])
    Wsq <- diag(sqrt(exp(fhat[,i])))
    L <- t(chol(diag(nrow=t)+Wsq%*%Kx%*%Wsq))
    three <- as.numeric(sum(log(diag(L))))
    sum <- sum+(one-two-three)
  }
  prior <- 0
  for (j in 1:q){
    prior <- prior +t(X[,j])%*%Ktinv%*%X[,j]
  }
  
  return((-1)*(sum-0.5*prior))
}
