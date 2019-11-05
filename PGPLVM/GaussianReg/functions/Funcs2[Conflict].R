#Basic vanilla GPLVM
#Optimise wrt x

#Loglik
loglikchol <- function(par){
  dim <- dim(y)[2]
  #Calculate covariance mat
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  x <- par[-(1:3)]
  K <- fillup(x,x,kern,par[1:3])
  #Find cholesky decomp
  L <- t(chol(K+par[3]^2*diag(n)))
  result <- 0
  
  for(j in 1:dim){
    alpha <- solve(t(L),solve(L,y[,j]))
    result <- result + (-0.5*t(y[,j])%*%alpha-sum(log(diag((L)))))
  }
  return((-1)*(result-0.5*t(x)%*%x))
}

difkern2 <- function(index,x,l,a){
  #n <- dim(x)[1]
  output <- matrix(0,nrow=n,ncol=n)
  vec <- x-x[index]
  applied <- a^2*exp(-0.5*vec^2/l^2)
  output[index,] <- -vec*applied
  output[,index] <- vec*applied
  output[index,index] <- 0
  return(output)
}
#Gradient
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

gradientChol2.0 <- function(par){
  len <- length(par)
  gvec <- numeric(len)
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  x <- par[-(1:3)]
  #Find expression for K inverse
  K <- fillup(x,x,kern,par[1:3])
  L <- t(chol(K+sigma^2*diag(n)))
  Linv <- solve(L)
  Kinv <- t(Linv)%*%Linv
  result <- 0
  for (i in 1:len){
    if (i == 1){
      dKi <- fillup(x,x,dKl,par[1:3])
    } else if (i == 2){
      dKi <- fillup(x,x,dKa,par[1:3])
    } else if (i == 3){
      dKi <- dKsigma(sigma)
    } else{
      dKi <- difkern2((i-3),x,l,a)
    }
    for (j in 1:dim){
      alpha <- solve(t(L),solve(L,y[,j]))
      result <- result +0.5*t(alpha)%*%dKi%*%alpha-0.5*sum(diag(Kinv%*%dKi))
    }
    gvec[i] <- result
    result <- 0
  }
  #Fix p(x) term
  gvec[-(1:3)] <- gvec[-(1:3)]-x
  gvec[1] <- (-1)*gvec[1]
  return(-gvec)
}


gradientchol <- function(par){
  dim <- dim(y)[2]
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  x <- par[-(1:3)]
  len <- length(par)
  K <- fillup(x,x,kern,par[1:3])
  t <- length(x)
  L <- t(chol(K+sigma^2*diag(n)))
  gvec <- numeric(len)
  result <- 0
  for (i in 1:len){
    dKx <- difkern(i,x,K)
    for (j in 1:dim){
      alpha <- solve(t(L),solve(L,y[,j]))
      Lint <- solve(L)
      Kinv <- t(Lint)%*%Lint
      result <- result +0.5*t(alpha)%*%dKx%*%alpha-0.5*sum(diag(Kinv%*%dKx))
    }
    gvec[i] <- result
    result <- 0
  }
  return(-gvec)
}
predictor <- function(xstar,x,dimm){
  dim <- 2
  #Calculate covariance mat
  K <- fillup(x,x,kern,l)
  #Find cholesky decomp
  L <- t(chol(K+sigma^2*diag(n)))
  
  kstar <- fillup(x,xstar,kern,l)
  alpha <- solve(t(L),solve(L,y[,dimm]))
  return(t(kstar)%*%alpha)
}
fillup <- function(x1,x2,kernel,theta){
  n1 <- length(x1)
  n2 <- length(x2)
  matr2 <- matrix(NA,nrow=n1,ncol=n2)
  for (i in 1:n1){
    for (j in 1:n2){
      matr2[i,j] <- kernel(x1[i],x2[j],theta)
    }
  }
  return(matr2)
}

kern <- function(x1,x2,theta){
  return(theta[2]^2*exp(-(x1-x2)^2/(2*theta[1]^2)))
}

dKl <- function(x1,x2,theta){
  return(theta[2]^2*exp(-(x1-x2)^2/(2*theta[1]^2))*(-(x1-x2)^2/(theta[1]^3)))
}
dKa <- function(x1,x2,theta){
  return(2*theta[2]*exp(-(x1-x2)^2/(2*theta[1]^2)))
}
dKsigma <- function(sigma){
  return(diag(n)*sigma*2)
}

