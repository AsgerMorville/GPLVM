Kx_obj <- function(x,model){
  x <- as.vector(x)
  x <- matrix(x,ncol=model$q,nrow=model$t)
  Kx <- fillup(x,ARDKernel,model$thetaf)+diag(model$t)*10e-05
  U <- chol(Kx)
  Uinv <- backsolve(r=U,x=diag(ncol(U)))
  Kxinv <- Uinv%*%t(Uinv)
  
  return(list(Kx=Kx,Kxinv=Kxinv))
}
