#Generate data

#Generate data function. Outputs a model-object
library(MASS)
generateData <- function(t,p,q){
  set.seed(1)
  model <- list()
  model$t <- t
  model$q <- q
  model$p <- p
  model$tmat <- matrix(1:t,ncol=1)
  model$thetax <- list(r=1,l=1)
  model$thetaf <- list(rho=1,delta=1)
  #model$beta <- 100
  
  tvec <- seq(from=1,to=t,length.out = t)
  Kt <- fillup(as.matrix(tvec),matern,model$thetax)
  model$x <- matrix(NA,nrow=t,ncol=q)
  for (l in 1:q){
    model$x[,l] <- mvrnorm(n=1,mu=rep(0,t),Sigma=Kt)
  }
  model$x <- scale(model$x,scale=F)
  Kx <- fillup(as.matrix(model$x),ARDKernel,model$thetaf)
  
  model$f <- matrix(NA,nrow=t,ncol=p)
  model$y <- matrix(NA,nrow=t,ncol=p)
  for (l in 1:p){
    model$f[,l] <- mvrnorm(n=1,mu=rep(0,t),Sigma=Kx)
    model$f <- scale(model$f,scale=F)
    model$y[,l] <- rpois(model$t,exp(model$f[,l]))
  }
  model$Ktinv <- solve(Kt)
  model$Kt <- Kt
  return(model)
}


