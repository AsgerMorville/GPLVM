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


x <- matrix(runif(200),nrow=100,ncol=2)

fillup2(x,kern,c(1,1,0.1))

#rfillup(x,x,kern,c(1,1,0.1))
