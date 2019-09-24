dPsi1dmu <- function(nind,qind,mod){
  n <- dim(mod$psi1)[1]
  m <- dim(mod$psi1)[2]
  muqn <- mod$mu[nind,qind]
  wj <- mod$thetaf$lengthscales[qind]
  S <- mod$Lambda[nind,qind]
  end <- matrix(0,nrow=n,ncol=m)
  xu <- mod$xu
  factor <- -(wj*(muqn-xu[,qind]))/(wj*S+1)
  asdf <- 123
  end[nind,] <- factor
  
  return(mod$psi1*end)
}

