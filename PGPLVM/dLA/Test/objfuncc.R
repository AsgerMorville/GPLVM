
objfunc <- function(vec){
  xmat <- matrix(vec,ncol=model$q,nrow=model$t)
  Kx_object <- Kx_obj(vec)
  Kxinv <- Kx_object$Kxinv
  Kx <- Kx_object$Kx
  
  #Skinv <- solve(model$S)
  #A <- Skinv+Kxinv
  #Ainv <- solve(A)
  
  
  sum <- 0
  for (i in 1:model$p){
    A <- Skinv_arr[,,i]+Kxinv
    Ainv <- solve(A)
    
    fhat <- as.vector(Ainv%*%model$A_arr[,,i]%*%model$fhatk[,i])
    first <- as.numeric(t(model$y[,i])%*%(fhat))-sum(exp(fhat))
    second <- 0.5*t(fhat)%*%Kxinv%*%fhat
    
    #Cholsum <- chol(diag(model$t)+Kx%*%Skinv_arr[,,i])
    
    #third <- sum(log(diag(Cholsum+diag(model$t)*10)))
    third <- 0.5*determinant(diag(model$t)+Kx%*%Skinv_arr[,,i],logarithm = T)$modulus
    sum <- sum+first-second-third
  }
  #Term involving Kt (prior)
  sumprior <- 0
  for (j in 1:model$q){
    sumprior <- sumprior + t(xmat[,j])%*%model$Ktinv%*%xmat[,j]
  }
  sumprior <- -0.5*sumprior
  return(as.numeric(-sum-sumprior))
}
