#Main dLa driver

dLA <- function(dataStart){
  #Input data y, start x
  #output: optimised x-positions. same dimension as startx matrix t times q
  y <- dataStart$y
  x <- dataStart$x0
  
  model <- list()
  model$t <- dim(x)[1]
  model$p <- dim(y)[2]
  model$q <- dim(x)[2]
  model$thetaf <- dataStart$thetaf
  model$y <- y
  maxIter <- dataStart$maxIter
  
  #Ktinv calculation
  tvec <- 1:model$t
  Kt <- fillup(as.matrix(tvec),matern,dataStart$thetax)
  model$Ktinv <- solve(Kt)
  output <- matrix(NA,ncol=model$q,nrow=model$t)
  for (a in 1:dataStart$maxiter){
    model$x <- x
    #Precomputations. We need to determine the expression to be maximized in step 5.
    #this will depend on 
    kxObj <- Kx_obj(x,model)
    model$Kxinv <- kxObj$Kxinv
    model$Kx <- kxObj$Kx
    
    #1. compute posterior mode given X. This will be a t times p matrix
    model$fhatk <- newton(y,model$Kxinv,100,model)
    
    #Also compute the A's
    A_arr <- array(NA,dim=c(model$t,model$t,model$p))
    for (i in 1:model$p){
      A_arr[,,i] <- diag(exp(model$fhatk[,i]))+kxObj$Kxinv
    }
    model$A_arr <- A_arr
    
    #2. Derive mk, sk
    Sk_arr <- array(NA,dim=c(model$t,model$t,model$p))
    Skinv_arr <- array(NA,dim=c(model$t,model$t,model$p))
    for (i in 1:model$p){
      Skinv_arr[,,i] <- A_arr[,,i]-kxObj$Kxinv
      Sk_arr[,,i] <- solve(Skinv_arr[,,i])
    }
    model$Skinv_arr <- Skinv_arr
    model$Sk_arr <- Sk_arr
    #Sk <- solve(A+Kxinv)
    #mk <- Sk%*%A%*%fhat
    
    #3. Find function to be updated and the gradient
    x <- xOptim(model)
    if (sum(abs(model$x-x))<0.1){break}
    #model$x <- x
  }
  output <- x
  return(output)
}


