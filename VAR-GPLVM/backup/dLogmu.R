#Gradient of log likelihood wrt. variational means
dLogmu <- function(psi,Kxx,Kuu,Kuuinv,y,sigma,p){
  #First calculate dFhat/dMuj for each j in q
  psi0 <- psi$psi0
  psi1 <- psi$psi1
  psi2 <- psi$psi2
  YYt <- y%*%t(y)
  A <- sigma^2*Kuu+psi2
  Ainv <- solve(A)
  psi1trace <- YYt%*%psi1%*%Ainv
  psi2trace <- p*Kuuinv-sigma^2*p*Ainv-Ainv%*%t(psi1)%*%YYt%*%psi1%*%Ainv
  
}