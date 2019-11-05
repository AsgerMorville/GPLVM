newton <- function(y,Kx_inv,maxiter,model){
  f <- matrix(0,ncol=model$p,nrow=model$t)
  for (u in 1:model$p){
    for(j in 1:maxiter){
      gradvec <- gradientNewton(f[,u],y[,u],Kx_inv)
      checker <- norm(gradvec,type="2")
      if (checker<0.001){
        break
      }
      invhMat <- solve(hessianNewton(f[,u],Kx_inv))
      f[,u] <- as.numeric(f[,u])-as.numeric(0.1*invhMat%*%gradvec)
    }
  }
  return(f)
}