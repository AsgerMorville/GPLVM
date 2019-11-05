#wrpaer

#Let gradient and loglikelihoods take vector as argument, and output is vector 

vecLogLik <- function(x){
  model$X <- matrix(x,ncol=model$q,byrow=F)
  #
  return(fullLogLik(model))
}

vecGrad <- function(x){
  #model <<- model
  model$X <- matrix(x,ncol=model$q,byrow=F)
  return(as.vector(fullgrad2(model)))
  
}
