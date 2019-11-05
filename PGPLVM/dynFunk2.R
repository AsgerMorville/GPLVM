#hardcode
theta <- c(1,1,0.1)

logLik <- function(x){
  t <- dim(x)[1]
  K <- fillup2(x,kern,theta)
  obj <- chol_inverse(K+sigma^2*diag(t))
  Kinv <- obj$Kinv
  L <- obj$L
  YYt <- (y)%*%t(y)
  return((-1)*(-dim*sum(log(diag((L))))-0.5*sum(diag(Kinv%*%YYt))-0.5*t(x)%*%Ktinv%*%x))
}
gradLogLik <- function(X){
  t <- dim(X)[1]
  q <- dim(X)[2]
  K <- fillup2(X,kern,theta)
  obj <- chol_inverse(K+sigma^2*diag(t))
  Kinv <- obj$Kinv
  L <- obj$L
  YYt <- y%*%t(y)
  dl <- Kinv%*%YYt%*%Kinv-dim*Kinv
  gvec <- matrix(NA,nrow=t,ncol=q)
  for (i in 1:t){
    for (j in 1:q){
      dKi <- difkernq(i,j,X,K,l)
      #dKi <- difkern(i,as.vector(X),K)
      result <- sum(diag(dl%*%dKi))
      gvec[i,j] <- result
      result <- 0
    }
  }
  return(-(gvec-t(t(X)%*%Ktinv)))
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
  q<-1
  n <- dim(x1)[1]
  mat <- array(NA,dim=c(n,n,q))
  for (i in 1:q){
    mat[,,i] <- (matrix(x1[,i],t,t,byrow=F)-matrix(x1[,i],t,t,byrow=T))^2
  }
  mat <- apply(mat,c(1,2),sum)
  fin <- apply(mat,c(1,2),kern,theta)
  return(fin)
}
fillup <- function(x1,x2,kern,l){
  n1 <- length(x1)
  n2 <- length(x2)
  matr2 <- matrix(NA,nrow=n1,ncol=n2)
  for (i in 1:n1){
    for (j in 1:n2){
      matr2[i,j] <- kern(x1[i],x2[j],l)
    }
  }
  return(matr2)
}
kern1 <- function(x1,x2,l){
  return(exp(-(x1-x2)^2/(2*l^2)))
}
kern <- function(r,theta){
  return(theta[2]^2*exp(-r/(2*theta[1]^2)))
}
chol_inverse <- function(A){
  L <- t(chol(A))
  Linv <- solve(L)
  return(list(Kinv=t(Linv)%*%Linv,L=L))
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
    if (noiter >250){break}
    
  }
  return(x1)
}

