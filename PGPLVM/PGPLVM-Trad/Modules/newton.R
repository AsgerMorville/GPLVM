newton <- function(y_i,X,graddy,Hessian,maxiter,Kx_inv){
  t <- dim(X)[1]
  f <-numeric(t)
  for(j in 1:maxiter){
    grad <-as.numeric(graddy(f,y_i,Kx_inv))
    checker <- norm(grad,type="2")
    if (checker<0.001){
      break
    }
    invH <- solve(Hessian(f,Kx_inv))
    f <- as.numeric(f)-as.numeric(0.1*invH%*%grad)
    
  }
  return(f)
}