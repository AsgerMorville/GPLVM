#vectoList

vecLogLik <- function(vecto){
  modeltest$hyper$alpha <- vecto[1]
  modeltest$hyper$gamma <- vecto[2]
  modeltest$hyper$beta <- vecto[3]
  
  #q <- dim(modeltest$Xny)[2]
  modeltest$x <- matrix(vecto[-c(1:3)],ncol=modeltest$q)
  
  return(logLik(modeltest))
}



vecGrad <- function(vector){
  modeltest$hyper$alpha <- vector[1]
  modeltest$hyper$gamma <- vector[2]
  modeltest$hyper$beta <- vector[3]
  #q <- dim(modeltest$Xny)
  modeltest$x <- matrix(vector[-c(1:3)],ncol=modeltest$q)
  obj <- gradLogLikv2(modeltest)
  #return(modeltest)
  return(c(obj$hyper,obj$x))
}
