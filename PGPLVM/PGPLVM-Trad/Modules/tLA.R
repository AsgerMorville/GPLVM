tLA <- function(dataArg){
  model <- list()
  model$thetaf <- dataArg$thetaf
  model$thetax <- dataArg$thetax
  
  model$x <- dataArg$x
  model$y <- dataArg$y
  model$q <- dim(model$x)[2]
  model$t <- dim(model$x)[1]
  model$p <- dim(model$y)[1]
  
  return(scg(vecLogLik,vecGrad,c(as.vector(model$x)),1))
}
