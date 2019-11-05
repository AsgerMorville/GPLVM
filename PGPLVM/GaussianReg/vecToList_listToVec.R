vecToList <- function(vec,t,q,p){
  #This function takes a vector as input, and treats the first part as 
  #1. variational means
  #2. Lambda values, bycol
  #Define varmeans
  mu <- matrix(vec[1:n*q],nrow=t,ncol=q,byrow=F)
  #Define lambda matrix
  Lambda <- matrix(vec[(n*q+1):(2*n*q)],nrow=t,ncol=q)
  return(list(mu=mu,Lambda=Lambda))
}
listToVec <- function(list){
  return(unlist(list))
}