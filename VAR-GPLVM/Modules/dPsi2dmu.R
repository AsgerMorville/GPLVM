dPsi2dmu <- function(nind,qind,mod){
  psi2 <- mod$psi2
  thetaf <- mod$thetaf
  mu <- mod$mu
  l <- mod$Lambda
  xu <- mod$xu
  w <- thetaf$lengthscales
  
  m <- dim(psi2)[1]
  muqn <- mu[nind,qind]
  wj <- as.numeric(w[qind])
  S <- as.numeric(l[nind,qind])
  
  #Find the nind'th psi2 matrix
  psi2n <- matrix(NA,ncol=m,nrow=m)
  Si <- as.vector(l[nind,])
  for (k in 1:m){
    for (kp in 1:m){
      xuk <- as.vector(xu[k,])
      xukp <- as.vector(xu[kp,])
      w <- as.vector(w)
      mui <- as.vector(mu[nind,])
      
      xbar <- (xuk+xukp)/2
      top1 <- w*(xuk-xukp)^2
      top2 <- w*(mui-xbar)^2
      bottom2 <- 2*w*Si+1
      Toppen <- exp(-top1/4-top2/bottom2)
      psi2n[k,kp] <- thetaf$sigmaf^4*prod(Toppen/sqrt(bottom2)) 
    }
  }
  
  
  xbarmat <- outer(xu[,qind],xu[,qind],'+')/2
  end2 <- -2*wj*(muqn-xbarmat)/(2*wj*S+1)
  return(psi2n*end2)
  #return(xbarmat)
}
