gradLogLikv2 <- function(modeltest){
  y <- modeltest$y
  X <- modeltest$x
  dim <- dim(y)[2]
  n <- dim(X)[1]
  q <- dim(X)[2]
  t <- n
  Ktinv <- modeltest$Ktinv
  
  a <- modeltest$hyper$alpha
  sigma <- 1/sqrt(modeltest$hyper$beta)
  thetaf <- list()
  thetaf$alpha <-modeltest$hyper$alpha
  thetaf$gamma <- modeltest$hyper$gamma
  #par <- c(modeltest$hyper$gamma,modeltest$hyper$alpha,1/sqrt(modeltest$hyper$beta))
  #l <- par[1]
  K <- fillup2(X,RBF_vec,thetaf)
  obj <- chol_inverse(K+sigma^2*diag(n))
  Kinv <- obj$Kinv
  L <- obj$L
  YYt <- y%*%t(y)
  dl <- 0.5*(Kinv%*%YYt%*%Kinv-dim*Kinv)
  gvec <- matrix(NA,nrow=n,ncol=q)
  for (i in 1:n){
    for (j in 1:q){
      dKi <- difkernq(i,j,X,K,thetaf$gamma)
      #dKi <- difkern(i,as.vector(X),K)
      result <- sum(diag(dl%*%dKi))
      gvec[i,j] <- result
    }
  }
  dKgamma <- fillup2(X,dKg,thetaf)
  dKalpha <-  fillup2(X,dKa,thetaf)
  dKsig <- dKsigma(modeltest$hyper$beta)
  dgamma <- sum(diag(dl%*%dKgamma))
  dalpha <- sum(diag(dl%*%dKalpha))
  dSig <- sum(diag(dl%*%dKsig))
  hyper <- c(-dalpha,-dgamma,-dSig)
  xGrad <- -(gvec-t(t(X)%*%Ktinv))
  return(list(x=xGrad,hyper=hyper))
}
