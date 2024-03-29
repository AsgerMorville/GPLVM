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
