dPsi2dxu <- function(kpp,j,mod){
  #First, find the unsummed psi2-array
  sigmaf <- thetaf$sigmaf
  
  psi2temp <- array(NA,dim=c(m,m,n))
  Lambda <-mod$Lambda
  xu <- mod$xu
  w <- mod$thetaf$lengthscales
  mu <- mod$mu
  m <- dim(mod$xu)[1]
  t <- dim(mod$mu)[1]
  for (i in 1:t){
    Si <- as.vector(Lambda[i,])
    #Si <- as.vector(S[i,i,])
    for (k in 1:m){
      for (kp in 1:m){
        xuk <- as.vector(xu[k,])
        xukp <- as.vector(xu[kp,])
        w <- as.vector(w)
        mui <- as.vector(mu[i,])
        
        xbar <- (xuk+xukp)/2
        top1 <- w*(xuk-xukp)^2
        top2 <- w*(mui-xbar)^2
        bottom2 <- 2*w*Si+1
        Toppen <- exp(-top1/4-top2/bottom2)
        psi2temp[k,kp,i] <- sigmaf^4*prod(Toppen/sqrt(bottom2)) 
      }
    }
  }
  wj <- w[j]
  xukj <- xu[kpp,j]
  end <- array(0,dim=c(m,m,t))
  #For each i, multiply by the factor(derivative)
  for (i in 1:t){
    Si <- Lambda[i,j]
    temp <- matrix(0,ncol=m,nrow=m)
    temp[kpp,] <- as.vector(-wj*(xukj-xu[,j])/2+wj*(2*mu[i,j]-xukj-xu[,j])/(2*wj*Si+1))
    temp[,kpp] <- as.vector(wj*(xu[,j]-xukj)/2+wj*(2*mu[i,j]-xukj-xu[,j])/(2*wj*Si+1))
    temp[kpp,kpp] <- as.numeric(wj*(2*mu[i,j]-xukj)/(2*wj*Si+1))
    end[,,i] <- psi2temp[,,i]*temp
  }
  
  return(apply(end,c(1,2),sum))
}
