#Log likelihood
#logLik <- function(par,X){
#  dim <- dim(y)[2]
#  l <- par[1]
#  a <- par[2]
#  q <- dim(X)[2]
#  sigma <- par[3]
#  K <- fillup2(X,RBF_vec,par[1:3])
#  obj <- chol_inverse(K+sigma^2*diag(t))
#  Kinv <- obj$Kinv
#  L <- obj$L
#  YYt <- y%*%t(y)
#  prior <- 0 
#  for (i in 1:q){
#    prior <- prior + t(x[,i])%*%Ktinv%*%x[,i]
#  }
#  prior <- -0.5*prior
#  return((-1)*(-dim*sum(log(diag((L))))-0.5*sum(diag(Kinv%*%YYt))-0.5*t(X)%*%Ktinv%*%X))
#}

#Log likelihood
logLik <- function(par,X){
  dim <- dim(y)[2]
  l <- par[1]
  a <- par[2]
  q <- dim(X)[2]
  sigma <- par[3]
  K <- fillup2(X,RBF_vec,par[1:3])
  obj <- chol_inverse(K+sigma^2*diag(t))
  Kinv <- obj$Kinv
  L <- obj$L
  YYt <- y%*%t(y)
  prior <- 0 
  for (i in 1:q){
    prior <- prior + t(X[,i])%*%Ktinv%*%X[,i]
  }
  prior <- -0.5*prior
  return((-1)*(-dim*sum(log(diag((L))))-0.5*sum(diag(Kinv%*%YYt))+prior))
}
