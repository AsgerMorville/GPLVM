# Functions for P-GPLVM

#Implementation of Pillow

optimiser <- function(y,xstart,maxiter){
  
  N <- dim(y)[1]
  t <- dim(y)[2]
  X <- xstart
  #1. a) Computer posterior mode of f_i for i in 1:N

  #Compute inverse
  Kx <- fillup(X,X,kern,1)
  L <- t(chol(Kx+diag(t)*10e-8))
  #Linv <- backsolve(L, diag(dim(L)[1]))
  Linv <- solve(L)
  Kx_inv <- t(Linv)%*%Linv
  
  fhat <- matrix(NA,nrow=N,ncol=t)
  for (i in 1:N){
    fhat[i,] <- newton(y[i,],X,gradient2,Hessian,maxiter,Kx_inv)
    print(i)
    }
  
  #b) Find Ak
  A <- array(NA,dim=c(t,t,N))
  for (i in 1:N){
    A[,,i] <- diag(exp(fhat[i,]))-Kx_inv
  }    
  
  #2. Derive m_k and S_k
  S <- array(NA,dim=c(t,t,N))
  m <- matrix(NA,nrow=t,ncol=N)
  for (i in 1:N){
    S[,,i] <- solve(A[,,i]-Kx_inv)
    m[,i] <- S[,,i]%*%A[,,i]%*%fhat[i,]
  }
  
  ## 5.  finding posterior X
  #Function to be minimised
  Xfunc <- function(x){
    Kx <- fillup(x,x,kern,1)
    L <- t(chol(Kx+diag(t)*10e-8))
    #Linv <- backsolve(L, diag(dim(L)[1]))
    Linv <- solve(L)
    Kx_inv <- t(Linv)%*%Linv
    result <- 0
    for (i in 1:N){
      #1st term
      one <- sum(y[i,]*fhat[i,])-sum(exp(fhat[i,]))-sum(log(factorial(y[i,])))
      one <- as.numeric(one)
      #2nd term
      two <- 0.5*t(fhat[i,])%*%Kx_inv%*%(fhat[i,])
      two <- as.numeric(two)
      #3rd
      Wsq <- diag(sqrt(exp(fhat[i,])))
      L <- t(chol(diag(nrow=t)+Wsq%*%Kx%*%Wsq))
      three <- as.numeric(sum(log(diag(L))))
      #eig <- eigen(diag(nrow=t)+Kx%*%W)
      #three <- 0.5*sum(log(eig$values))
      result <- result + (one-two-three)
    }
    return(result)
  }
  #return(Xfunc(xstart))
  #W <- -diag(exp(fhat[1,]))
  #L <- eigen(diag(nrow=t)+Kx%*%W)
  #return(L)
  return(fhat)
}




fullGrad <- function(tvec,y,xstart,maxiter){
  
  N <- dim(y)[1]
  t <- dim(y)[2]
  X <- xstart
  Kt <- fillup(tvec,tvec,matern,20)
  Ktinv <- solve(Kt)
  #1. a) Computer posterior mode of f_i for i in 1:N
  
  #Compute inverse
  Kx <- fillup(X,X,kern,1)
  L <- t(chol(Kx+diag(t)*10e-8))
  #Linv <- backsolve(L, diag(dim(L)[1]))
  Linv <- solve(L)
  Kx_inv <- t(Linv)%*%Linv
  
  fhat <- matrix(NA,nrow=N,ncol=t)
  for (i in 1:N){
    fhat[i,] <- newton(y[i,],X,gradientNewton,Hessian,maxiter,Kx_inv)
  }
  return((-1)*gradient8919(X,y,fhat,Kx,Kx_inv,Ktinv))
}

gradient8919 <- function(x,y,fhat,Kx,Kxinv,Ktinv){
  #Returns the gradient. Output format: matrix. Input: x = matrix.
  
  N <- dim(y)[1]
  t <- dim(y)[2]
  output <- matrix(0,nrow=t,ncol=1)
  
  for (i in 1:t){
    #Loop over the X_i coordinates.
    dKx <- difkernq(i,1,x,Kx,1)
    for (j in 1:N){
      first <- 0.5*fhat[j,]%*%Kxinv%*%dKx%*%Kxinv%*%fhat[j,]
      Wj <- diag(exp(fhat[j,]))
      inv <- solve(diag(t)+Kx%*%Wj)
      prod <- dKx%*%Wj
      second <- 0.5*sum(diag(inv%*%prod))
      output[i] <- output[i]+ first-second
    }
  }
  
  return(output-Ktinv%*%x)
}

gradientaLa <- function(x){
  sum <- 0
  for (i in 1:N){
    dKx <- difkernq(i,1,X,Kx,1)
    first <- 0.5*fhat[i,]%*%Kxinv%*%dKx%*%Kxinv%*%fhat[i,]
    inv <- solve(diag(t)+K%*%Wi)
    prod <- dKx%*%Wi
    second <- -0.5*sum(diag(inv%*%prod))
    sum <- sum+first+second
  }
  return(sum-t(t(X)%*%Ktinv))
  
}
matern <- function(x1,x2,l){
  return(exp(-abs(x1-x2)/l))
}



gradientNewton <- function(f,y,Kx_inv){
  return(y-exp(f)-Kx_inv%*%f)
}
Hessian <- function(f,Kx_inv){
  return(diag(-exp(f))-Kx_inv)
}

kern <- function(r,l){
  return(exp(-(r)^2/(2*l^2)))
}

kern1 <- function(x1,x2,l){
  return(exp(-(x1-x2)^2/(2*l^2)))
}
difkernq <- function(index,j,x,K,l){
  n <- dim(x)[1]
  x <- x[,j]
  output <- matrix(0,nrow=n,ncol=n)
  vec <- x[index]-x
  output[index,] <- K[index,]
  output[,index] <- K[,index]
  output[index,] <- output[index,]*(x-x[index])/l^2
  output[,index] <- output[,index]*(x-x[index])/l^2
  
  output[index,index] <- 0
  return(output)
}
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
  return(x1)
}
