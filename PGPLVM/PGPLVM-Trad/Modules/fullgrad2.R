fullgrad2 <- function(model){
  #model contains tvec,y,xstart,maxiter.
  #First, fhat is found. Then we fix fhat and calculate derivative.
  y <- model$y
  tvec <- model$tvec
  X <- model$X
  maxiter <- model$maxiter
  gradopt <- model$gradopt
  
  N <- dim(y)[2]
  t <- dim(y)[1]
  
  Ktinv <- model$Ktinv
  #1. a) Computer posterior mode of f_i for i in 1:N
  
  #Compute inverse
  Kx <- fillup(X,ARDKernel,1)
  L <- t(chol(Kx+diag(t)*10e-8))
  #Linv <- backsolve(L, diag(dim(L)[1]))
  Linv <- solve(L)
  Kx_inv <- t(Linv)%*%Linv
  
  #Compute fhat values
  fhat <- matrix(NA,nrow=t,ncol=N)
  for (i in 1:N){
    fhat[,i] <- newton(y[,i],X,gradientNewton,Hessian,maxiter,Kx_inv)
  }
  
  #Calculate gradient
  mod2 <- list()
  mod2$X <- X
  mod2$y <- y
  mod2$fhat <- fhat
  mod2$Kx <- Kx
  mod2$Kx_inv <- Kx_inv
  mod2$Ktinv <- Ktinv
  if (gradopt == "aLa"){
    finGrad <- gradient2(mod2)
  }
  if (gradopt == "tLa"){
    finGrad <- gradientTLA(mod2)
  }
  
  return((-1)*finGrad)
}
