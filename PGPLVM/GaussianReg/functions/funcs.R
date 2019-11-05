##Functions

grad5_9 <- function(par){
  K_y <- fillup(x,x,kern,par[1])+par[2]^2*diag(n)
  invx <- solve(K_y)
  yK <- t(y)%*%invx
  Ky <- invx%*%y
  pK <- fillup(x,x,difKern2,par[1])
  First <- 0.5*yK%*%pK%*%Ky
  second <- 0
  for(i in 1:n){
    second <- second + invx[i,]%*%pK[,i]
  }
  second <- second*(-0.5)
  return((-1)*c(First+second,2*par[2]*0.5*yK%*%Ky-0.5*2*par[2]*sum(diag(invx))))
}
simpelGrad <- function(par){
    K_y <- fillup(x,x,kern,par)+0.01*diag(n)
    invx <- solve(K_y)
    yK <- t(y)%*%invx
    Ky <- invx%*%y
    pK <- fillup(x,x,difKern2,par)
    First <- 0.5*yK%*%pK%*%Ky
    second <- 0
    for(i in 1:n){
      second <- second + invx[i,]%*%pK[,i]
    }
    second <- second*(-0.5)
    return(-c(First+second))
  }

difKern2 <- function(x1,x2,l){
  return(exp(-(x1-x2)^2/(2*l^2))*(-(x1-x2)^2/(l^3)))
}

loglik <- function(par){
  Kernel <- fillup(x,x,kern,par[1])
  invx <- solve(Kernel+diag(nrow=n)*par[2]^2)
  return((-1)*(-0.5*t(y)%*%invx%*%y-0.5*log(abs(det(Kernel+diag(nrow=n)*par[2]^2)))))
}
simpelloglik <- function(par){
  Kernel <- fillup(x,x,kern,par)
  invx <- solve(Kernel+diag(nrow=n)*0.01)
  return(-(-0.5*t(y)%*%invx%*%y-0.5*log(abs(det(Kernel+diag(nrow=n)*0.01)))))
}

Cholpred <- function(par,test){
  K <- fillup(x,x,kern,par[1])
  L <- t(chol(K+par[2]^2*diag(n)))
  alpha <- solve(t(L),solve(L,y))
  Kstar <- fillup(test,x,kern,par[1])
  return(Kstar%*%alpha)
}

GAscent <- function(grad,init,gamma){
  a <- init
  n3 <- 500
  gsize <- numeric(n3)
  for (i in 1:n3){
    
    nab <- grad(a)
    gsize[i] <- norm(nab,type = "2")
    a <- a+gamma*nab
    
    #if (norm(grad(a),type="2")<0.01) {break}
    print(i)
  }
  return(list(a=a,gsizes=gsize))
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

kern <- function(x1,x2,l){
  return(exp(-(x1-x2)^2/(2*l^2)))
}

