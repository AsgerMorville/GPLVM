#Final GPLVM functions

logLik <- function(par,X){
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  K <- fillup2(X,kern,par[1:3])
  obj <- chol_inverse(K+sigma^2*diag(n))
  Kinv <- obj$Kinv
  L <- obj$L
  YYt <- y%*%t(y)
  return((-1)*(-dim*sum(log(diag((L))))-0.5*sum(diag(Kinv%*%YYt))))
}

gradLogLik <- function(par,X){
  n <- dim(X)[1]
  q <- dim(X)[2]
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  K <- fillup2(x,kern,par[1:3])
  obj <- chol_inverse(K+sigma^2*diag(n))
  Kinv <- obj$Kinv
  L <- obj$L
  YYt <- y%*%t(y)
  dl <- Kinv%*%YYt%*%Kinv-dim*Kinv
  gvec <- matrix(NA,nrow=n,ncol=q)
  for (i in 1:n){
    for (j in 1:q){
      dKi <- difkernq(i,j,X,K,l)
      result <- sum(diag(dl%*%dKi))
      gvec[i,j] <- result
      result <- 0
    }
  }
  return(-gvec)
}

chol_inverse <- function(A){
  L <- t(chol(A))
  Linv <- solve(L)
  return(list(Kinv=t(Linv)%*%Linv,L=L))
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
fillup2 <- function(x1,kern,theta){
  n <- dim(x1)[1]
  q <- dim(x1)[2]
  mat <- array(NA,dim=c(n,n,q))
  for (i in 1:q){
    mat[,,i] <- (matrix(x1[,i],n,n,byrow=F)-matrix(x1[,i],n,n,byrow=T))^2
  }
  mat <- apply(mat,c(1,2),sum)
  fin <- apply(mat,c(1,2),kern,theta)
  return(fin)
}

kern <- function(r,theta){
  return(theta[2]^2*exp(-r/(2*theta[1]^2)))
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
    x <- x1
    print(x[1:3])
    print(noiter)
    
  }
  return(x1)
}
