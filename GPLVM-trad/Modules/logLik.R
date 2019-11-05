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
#logLik <- function(par,X){
#  dim <- dim(y)[2]
#  l <- par[1]
#  a <- par[2]
#  q <- dim(X)[2]
#  sigma <- par[3]
#  K <- fillup2(X,RBF_vec,par[1:3])
# obj <- chol_inverse(K+sigma^2*diag(t))
#  Kinv <- obj$Kinv
#  L <- obj$L
#  YYt <- y%*%t(y)
#  prior <- 0 
#  for (i in 1:q){
#   prior <- prior + t(X[,i])%*%Ktinv%*%X[,i]
#  }
#  prior <- -0.5*prior
#  return((-1)*(-dim*sum(log(diag((L))))-0.5*sum(diag(Kinv%*%YYt))+prior))
#}
logLik <- function(model){
  thetaf <- list()
  thetaf$alpha <- model$hyper$alpha
  thetaf$gamma <- model$hyper$gamma
  #thetaf$beta <- 
  x <- model$x
  y <- model$y
  D <- dim(y)[2]
  gamma <- model$hyper$gamma
  alpha <- model$hyper$alpha
  Ktinv <- model$Ktinv
  YYt <- model$YYt
  q <- dim(x)[2]
  beta <- model$hyper$beta
  sigma <- 1/sqrt(model$hyper$beta)
  t <- dim(x)[1]
  
  #Precomputations
  K <- fillup2(x,RBF_vec,thetaf)
  
  obj <- chol_inverse(K+sigma^2*diag(t))
  Kinv <- obj$Kinv
  L <- obj$L
  #Computer prior term
  prior <- 0 
  for (i in 1:q){
    prior <- prior + t(x[,i])%*%Ktinv%*%x[,i]
  }
  prior <- -0.5*prior
  return((-1)*(-D*sum(log(diag((L))))-0.5*sum(diag(Kinv%*%YYt))+prior))
}
