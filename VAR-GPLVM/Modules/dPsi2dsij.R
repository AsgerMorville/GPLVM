dPsi2dsij <- function(i,j,mod){
  #Output is m by m
  m <- dim(mod$xu)[1]
  output <- matrix(NA,nrow=m,ncol=m)
  wj <- mod$thetaf$lengthscales[j]
  xu <- mod$xu
  mu <- mod$mu
  thetaf <- mod$thetaf
  for (k in 1:m){
    for (kp in 1:m){
      output[k,kp] <- mod$psi2[k,kp]*(2*(wj*(2*mod$mu[i,j]-mod$xu[k,j]-mod$xu[kp,j])/(2*(2*wj*mod$S_mat[i,j]+1)))^2-0.5*2*wj/(2*wj*mod$S_mat[i,j]+1))
    }
  }
  
  #Find the nind'th psi2 matrix
  psi2n <- matrix(NA,ncol=m,nrow=m)
  Si <- as.vector(mod$Lambda[i,])
  for (k in 1:m){
    for (kp in 1:m){
      xuk <- as.vector(xu[k,])
      xukp <- as.vector(xu[kp,])
      w <- wj
      mui <- as.vector(mu[i,])
      
      xbar <- (xuk+xukp)/2
      top1 <- w*(xuk-xukp)^2
      top2 <- w*(mui-xbar)^2
      bottom2 <- 2*w*Si+1
      Toppen <- exp(-top1/4-top2/bottom2)
      psi2n[k,kp] <- thetaf$sigmaf^4*prod(Toppen/sqrt(bottom2)) 
    }
  }
  return(psi2n*output)
}
