#Function to initialise the L-values

calibrateLambda <- function(n,q,tvec){
  vdcrange <- seq(from=0.05,to=3,by=1)
  initCovarMedian <- 0.18
  bestVdc <- 0.05
  bestDiff <- 1000
  
  
  for (vdc in vdcrange){
    #covars <- matrix(0.1+0.001*rnorm(n*q,mean=0,sd=1),nrow=n,ncol=q)
    #covars <- pmax(covars,0.05)
    covars <- vdc*matrix(1,nrow=n,ncol=q)
    Sj <- Sjcalculator(covars,tvec,n,q)
    
    med <- median(Sj)
    curdiff <- abs(med-initCovarMedian)
    if (curdiff < bestDiff){
      bestVdc <- vdc
      bestDiff <- curdiff
      bestMed <- med
    }
  }
  return(list(covars=bestVdc*matrix(1,nrow=n,ncol=q),vdc=bestVdc))
}

Sjcalculator <- function(L,tvec,n,q){
  Kx <- fillup(tvec,kern=matern,theta=modeltest$thetax)
  Kxinv <- solve(Kx)
  output <- array(0,dim=c(n,n,q))
  output2 <- matrix(NA,nrow=n,ncol=q)
  for (j in 1:q){
    output[,,j] <- solve(as.matrix(Kxinv+diag(as.vector(L[,j]))))
    output2[,j] <- pmax(diag(as.matrix(output[,,j])),0.05)
  }
  #Extract diagonals

  return(as.matrix(output2))
}
