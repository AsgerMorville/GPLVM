#dPsi1dwj

dPsi2dwj <- function(j,mod){
  psi1 <- mod$psi1
  m <- dim(mod$xu)[1]
  t <- dim(mod$y)[1]
  factor <- matrix(NA,ncol=m,nrow=t)
  mu <- mod$mu
  S_mat <- mod$S_mat
  xu <- mod$xu
  wj <- mod$thetaf$lengthscales[j]
  w <-  mod$thetaf$lengthscales
  sigmaf <-mod$thetaf$sigmaf
  
  output <- array(NA,dim=c(m,m,t))
  psi2temp <- array(NA,dim=c(m,m,t))
  #We create the m by m by n array of psi2-matrices. For each entry, we multiply with the diff.term
  for (i in 1:t){
    Si <- as.vector(S_mat[i,])
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
        
        #Find the term with which we have to multiply:
        term <- -(xu[k,j]-xu[kp,j])^2/4-((2*mu[i,j]-xu[k,j]-xu[kp,j])/(2*(2*wj*S_mat[i,j]+1)))^2-S_mat[i,j]/(2*wj*S_mat[i,j]+1)
        output[k,kp,i] <- psi2temp[k,kp,i]*term
      }
    }
  }
  return(apply(output,c(1,2),sum))
}
