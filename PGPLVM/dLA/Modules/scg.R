scg <- function(f,derf,x,tol){
  n <- length(x)
  deltastep <- 0
  lambda <- 1e-8
  lambdabar <- 0
  sigmac <- 1e-5
  sucess <- 1
  
  #Calculate initial gradient
  noiter <- 0
  
  pv <- -derf(x) 
  rv <- pv
  
  while(norm(rv,type="2")>tol){
    noiter <- noiter+1
    if(deltastep==0){
      df <- derf(x)
    } else {
      df <- -rv
    }
    print(norm(df,type="2"))
    deltastep <- 0
    if (sucess == 1){
      sigma <- sigmac/norm(pv,type="2")
      dfplus <- derf(x+sigma*pv)
      stilda <- (dfplus-df)/sigma
      delta <- as.numeric(t(pv)%*%stilda)
    }
    # Scale
    delta <- as.numeric(delta+(lambda-lambdabar)*norm(pv,type="2")^2)
    if (delta <= 0){
      lambdabar <- 2*(lambda-delta/norm(pv,type="2")^2)
      delta <- -delta+lambda*norm(pv,type="2")^2
      lambda <- lambdabar
    }
    # Step size
    mu <- as.numeric(t(pv)%*%rv)
    alpha <- mu/delta
    fv <- f(x)
    fvplus <- f(x+alpha*pv)
    delta1 <- as.numeric(2*delta*(fv-fvplus)/mu^2)
    rvold <- rv
    pvold <- pv
    if (delta1 >= 0){
      deltastep <- 1
      x1 <- x+alpha*pv
      rv <- -derf(x1)
      lambdabar <- 0
      sucess <- 1
      if (noiter%%n == 0){
        pv <- rv
      } else{
        rdiff <- rv-rvold
        beta <- as.numeric((t(rdiff)%*%rv)/as.numeric(t(rvold)%*%rvold))
        pv <- rv+beta*pvold
      }
      if (delta1 >= 0.75){
        lambda <- 0.25*lambda
      }
    } else{
      lambdabar<-lambda
      sucess <- 0
      x1 <- x+alpha*pv
    }
    if (delta1<0.25){
      lambda <- lambda+delta*(1-delta1)/norm(pvold,type="2")^2
    }
    if (norm(x-x1,type="2")<0.000001){break}
    x <- x1
    print(noiter)
    
    if (noiter >100){break}
    
  }
  return(x)
}
