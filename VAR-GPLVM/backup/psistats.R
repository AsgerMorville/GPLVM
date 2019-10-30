psistats <- function(mu,Lambda,xu,thetaf,n,m){
  psi0 <- n
  psi1 <- matrix(NA,nrow=n,ncol=m)
  psi2temp <- array(NA,dim=c(m,m,n))
  hypw <- thetaf$lengthscales
  sigmaf <- thetaf$sigmaf
  for (i in 1:n){
    Si <- as.vector(Lambda[i,])
    #Si <- as.vector(S[i,i,])
    for (k in 1:m){
      w <- as.vector(hypw)
      mui <- as.vector(mu[i,])
      xuk <- as.vector(xu[k,])
      top <- w*(mui-xuk)^2
      bot <- w*Si+1
      Fin <- prod(exp(-0.5*top/bot)/sqrt(bot))
      
      #top <- as.vector(hypw*as.vector(mu[i,]-xu[k,])^2)
      #bottom <- as.vector(2*(hypw*Si+1))
      #first <- exp(-top/bottom)
      
      #psi1[i,k] <- sigmaf^2*prod(first/sqrt(bottom))
      psi1[i,k] <- Fin
    }
  }
  
  for (i in 1:n){
    Si <- as.vector(Lambda[i,])
    #Si <- as.vector(S[i,i,])
    for (k in 1:m){
      for (kp in 1:m){
        xuk <- as.vector(xu[k,])
        xukp <- as.vector(xu[kp,])
        w <- as.vector(hypw)
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
  
  psi2 <- apply(psi2temp,c(1,2),sum)
  return(list(psi0=psi0,psi1=psi1,psi2=psi2))
  #return(list(psi0=psi0,psi1=psi1))
}

