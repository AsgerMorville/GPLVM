gradLogLik <- function(par,X){
  dim <- dim(y)[2]
  n <- dim(X)[1]
  q <- dim(X)[2]
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  K <- fillup2(X,RBF_vec,par[1:3])
  obj <- chol_inverse(K+sigma^2*diag(n))
  Kinv <- obj$Kinv
  L <- obj$L
  YYt <- y%*%t(y)
  dl <- Kinv%*%YYt%*%Kinv-dim*Kinv
  gvec <- matrix(NA,nrow=n,ncol=q)
  for (i in 1:n){
    for (j in 1:q){
      dKi <- difkernq(i,j,X,K,l)
      #dKi <- difkern(i,as.vector(X),K)
      result <- sum(diag(dl%*%dKi))
      gvec[i,j] <- result
      result <- 0
    }
  }
  dKlength <- fillup2(X,dKl,par)
  dKw <-  fillup2(X,dKa,par)
  dKsig <- dKsigma(sigma)
  dL <- sum(diag(dl%*%dKlength))
  dW <- sum(diag(dl%*%dKw))
  dSig <- sum(diag(dl%*%dKsig))
  hyper <- c(dL,dW,dSig)
  xGrad <- (gvec-t(t(X)%*%Ktinv))
  return(list(x=xGrad,hyper=hyper))
}
