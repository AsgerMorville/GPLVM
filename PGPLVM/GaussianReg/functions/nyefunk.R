#GPLVM functions
logNeil <- function(par){
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  x <- par[-(1:3)]
  K <- fillup2(x,kern,par[1:3])
  L <- t(chol(K+par[3]^2*diag(n)))
  Linv <- solve(L)
  Kinv <- t(Linv)%*%Linv
  yyt <- y%*%t(y)
  return((-1)*(-dim*sum(log(diag((L))))-0.5*sum(diag(Kinv%*%yyt))))
}

Ngrad <- function(par,X){
  if (is.vector(X)){
    n <- length(X)
    q <- 1
  }
  if (is.matrix(X)){
    n <- dim(X)[1]
    q <- dim(X)[2]
  }
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  K <- fillup2(x,kern,par[1:3])
  L <- t(chol(K+sigma^2*diag(n)))
  Linv <- solve(L)
  Kinv <- t(Linv)%*%Linv
  yyt <- y%*%t(y)
  dl <- Kinv%*%yyt%*%Kinv-dim*Kinv
  gvec <- numeric(n*q)
  if (q==1){
    for (i in 1:n){
        dKi <- difkern2(i,x,K,l)
        result <- sum(diag(dl%*%dKi))
        gvec[i] <- result
        result <- 0
    }
  }else{
    for (i in 1:n){
      for (j in 1:q){
        dKi <- difkernq(i,j,x,K,l)
        result <- sum(diag(dl%*%dKi))
        gvec[i] <- result
        result <- 0
      }
    }
  }
  
  return(gvec)
}
  
NeilGradient <- function(par){
  len <- length(par)
  gvec <- numeric(len)
  l <- par[1]
  a <- par[2]
  sigma <- par[3]
  x <- par[-(1:3)]
  #Find expression for K inverse
  K <- fillup2(x,kern,par[1:3])
  L <- t(chol(K+sigma^2*diag(n)))
  Linv <- solve(L)
  Kinv <- t(Linv)%*%Linv
  result <- 0
  for (i in 1:len){
    if (i == 1){
      dKi <- fillup2(x,dKl,par[1:3])
    } else if (i == 2){
      dKi <- fillup2(x,dKa,par[1:3])
    } else if (i == 3){
      dKi <- dKsigma(sigma)
    } else{
      dKi <- difkern2((i-3),x,K,l)
    }
    
    #for (j in 1:dim){
    #  alpha <- solve(t(L),solve(L,y[,j]))
    #  result <- result +0.5*t(alpha)%*%dKi%*%alpha-0.5*sum(diag(Kinv%*%dKi))
    #}
    if(i>3){
      yyt <- y%*%t(y)
      dl <- Kinv%*%yyt%*%Kinv-dim*Kinv
      result <- sum(diag(dl%*%dKi))
    }
    gvec[i] <- result
    result <- 0
  }
  #Fix p(x) term
  #gvec[-(1:3)] <- gvec[-(1:3)]-x
  gvec[1] <- (-1)*gvec[1]
  return(-gvec)
}

difkern2 <- function(index,x,K,l){
  #n <- dim(x)[1]
  n <- length(x)
  output <- matrix(0,nrow=n,ncol=n)
  vec <- x[index]-x
  #applied <- a^2*exp(-0.5*vec/l^2)
  output[index,] <- K[index,]
  output[,index] <- K[,index]
  output[index,] <- output[index,]*(x-x[index])
  output[,index] <- output[,index]*(x-x[index])
  
  output[index,index] <- 0
  return(output)
}
difkernq <- function(index,j,x,K,l){
  n <- dim(x)[1]
  #n <- length(x)
  x <- x[,j]
  output <- matrix(0,nrow=n,ncol=n)
  vec <- x[index]-x
  #applied <- a^2*exp(-0.5*vec/l^2)
  output[index,] <- K[index,]
  output[,index] <- K[,index]
  output[index,] <- output[index,]*(x-x[index])
  output[,index] <- output[,index]*(x-x[index])
  
  output[index,index] <- 0
  return(output)
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
#fillup <- function(x1,x2,kernel,theta){
#  n1 <- length(x1)
#  n2 <- length(x2)
#  matr2 <- matrix(NA,nrow=n1,ncol=n2)
#  for (i in 1:n1){
#    for (j in 1:n2){
#      matr2[i,j] <- kernel(x1[i],x2[j],theta)
#    }
#  }
#  return(matr2)
#}
fillup2 <- function(x1,kern,theta){
  if (is.vector(x1)){
    n <- length(x1)
    q <- 1
  }
  if (is.matrix(x1)){
    n <- dim(x1)[1]
    q <- dim(x1)[2]
  }
  mat <- array(NA,dim=c(n,n,q))
  if (q==1){
    mat[,,1] <- (matrix(x1,n,n,byrow=F)-matrix(x1,n,n,byrow=T))^2
  }else{
    for (i in 1:q){
      mat[,,i] <- (matrix(x1[,i],n,n,byrow=F)-matrix(x1[,i],n,n,byrow=T))^2
      
    }
    mat <- apply(mat,c(1,2),function(x){return(sum(x))})
    
  } 
  theta <- theta
  fin <- apply(mat,c(1,2),kern,theta)
  return(fin)
}

kern <- function(r,theta){
  return(theta[2]^2*exp(-r/(2*theta[1]^2)))
}
#kern <- function(x1,x2,theta){
#return(theta[2]^2*exp(-(norm(x1-x2,type="2")^2/(2*theta[1]^2))))
#}

dKl <- function(r,theta){
  return(theta[2]^2*exp(-r/(2*theta[1]^2))*(-r/(theta[1]^3)))
}
dKa <- function(r,theta){
  return(2*theta[2]*exp(-r/(2*theta[1]^2)))
}
dKsigma <- function(sigma){
  return(diag(n)*sigma*2)
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