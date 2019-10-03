vecToList <- function(vec,t,q,p,m){
  #This function takes a vector as input, and treats the first part as 
  #1. variational means
  #2. Lambda values, bycol
  #Define varmeans
  mubar <- matrix(vec[1:t*q],nrow=t,ncol=q,byrow=F)
  #Define lambda matrix
  Lambda <- matrix(vec[(t*q+1):(2*t*q)],nrow=t,ncol=q,byrow = F)
  #xu <- matrix(vec[(2*t*q+1):(2*t*q+m*q)],nrow=m,ncol=q,byrow=F)
  return(list(mubar=mubar,Lambda=Lambda))
  #,xu=xu))
  #return(list(mubar=mubar))
}
listToVec <- function(list){
  return(unlist(list))
}

vecLogLik <- function(vec){
  #Process vec so that it becomes on model-form
  #
  list <- vecToList(vec,t,q,p,m)

  modeltest$mubar <- list$mubar
  modeltest$Lambda <- list$Lambda
  #modeltest$xu <- list$xu
  return(-logLik(modeltest))
}
vecGradLogLik <- function(vec){
  #Process vec so that it becomes on model-form
  #
  list <- vecToList(vec,t,q,p,m)
  
  modeltest$mubar <- list$mubar
  modeltest$Lambda <- list$Lambda
  #modeltest$xu <- list$xu
  obj <- gradLoglik(modeltest)
  #Process obj so that it returns vec
  
  vec <- listToVec(obj)
  return(-vec)
}
