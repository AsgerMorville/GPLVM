#Loglik
loglikchol <- function(x){
  #Calculate covariance mat
  K <- fillup(x,x,kern,l)
  #Find cholesky decomp
  L <- t(chol(K+sigma^2*diag(n)))
  result <- 0
  for(j in 1:N){
    alpha <- solve(t(L),solve(L,yy[j,]))
    result <- result + (-0.5*t(yy[j,])%*%alpha-sum(log(diag((L)))))
  }
  return((-1)*(result-0.5*sum(t(x)%*%x)))
}


difkern <- function(index,x,K){
  t <- length(x)
  dKx <- matrix(0,ncol=t,nrow=t)
  dKx[index,] <- K[index,]
  dKx[,index] <- K[,index] 
  dKx[index,index] <- 0
  tempmat <- matrix(NA,nrow=t,ncol=t)
  for (k in 1:t){
    for (s in 1:t){
      tempmat[k,s] <- x[k]-x[s]
    }
  }
  tempmat <- tempmat/l^2
  dKx <- dKx*tempmat
  dKx[index,] <-  dKx[index,]*(-1)
  return(dKx)
}

gradientchol <- function(x){
  K <- fillup(x,x,kern,l)
  t <- length(x)
  L <- t(chol(K+sigma^2*diag(n)))
  gvec <- numeric(t)
  result <- 0
  for (i in 1:t){
    dKx <- difkern(i,x,K)
    for (j in 1:N){
      alpha <- solve(t(L),solve(L,yy[j,]))
      Lint <- solve(L)
      Kinv <- t(Lint)%*%Lint
      result <- result +0.5*t(alpha)%*%dKx%*%alpha-0.5*sum(diag(Kinv%*%dKx))
    }
    gvec[i] <- result
    result <- 0
    print(i)
  }
  return(-gvec)
}
