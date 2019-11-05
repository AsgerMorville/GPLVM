xOptim <- function(model){
  #function solving argmax_x p(y|x)p(x) (approximately)
  
  #Create obj. function
  objfunc <- function(vec){
    xmat <- matrix(vec,ncol=model$q,nrow=model$t)
    Kx_object <- Kx_obj(vec,model)
    Kxinv <- Kx_object$Kxinv
    Kx <- Kx_object$Kx
    sum <- 0
    for (i in 1:model$p){
      A <- model$Skinv_arr[,,i]+Kxinv
      Ainv <- solve(A)
      
      fhat <- as.vector(Ainv%*%model$A_arr[,,i]%*%model$fhatk[,i])
      first <- as.numeric(t(model$y[,i])%*%(fhat))-sum(exp(fhat))
      second <- 0.5*t(fhat)%*%Kxinv%*%fhat
      third <- 0.5*determinant(diag(model$t)+Kx%*%model$Skinv_arr[,,i],logarithm = T)$modulus
      sum <- sum+first-second-third
    }
    #Term involving Kt (prior)
    sumprior <- 0
    for (j in 1:model$q){
      sumprior <- sumprior + t(xmat[,j])%*%model$Ktinv%*%xmat[,j]
    }
    sumprior <- -0.5*sumprior
    return(as.numeric(-sum-sumprior))
  }
  
  #Create gradient
  gradObj <- function(vec){
    output <- matrix(0,ncol=model$q,nrow=model$t)
    xmat <- matrix(vec,ncol=model$q,nrow=model$t)
    Kx_object <- Kx_obj(vec,model)
    Kxinv <- Kx_object$Kxinv
    Kx <- Kx_object$Kx
    
    A_inv_x <- array(NA,dim=c(model$t,model$t,model$p))
    AKx_inv <- array(NA,dim=c(model$t,model$t,model$p))
    KxA_inv <- array(NA,dim=c(model$t,model$t,model$p))
    for (i in 1:model$p){
      A_inv_x[,,i] <- solve(as.matrix(model$Skinv_arr[,,i])+Kxinv)
      AKx_inv[,,i] <- A_inv_x[,,i]%*%Kxinv
      KxA_inv[,,i] <- Kxinv%*%A_inv_x[,,i]
    }
    for (l in 1:model$p){
      A <- model$Skinv_arr[,,l]+Kxinv
      Ainv <- solve(A)
      fhat <- as.vector(Ainv%*%model$A_arr[,,l]%*%model$fhatk[,l])
      inv <- solve(diag(model$t)+Kx%*%model$Skinv_arr[,,l])
      for (j in 1:model$q){
        for (i in 1:model$t){
          dKxdijmat <- dKxdxij(i,j,Kx,model$thetaf,xmat)
          dAinvij <- AKx_inv[,,l]%*%dKxdijmat%*%KxA_inv[,,l]
          dfhatl <- as.vector((dAinvij%*%model$A_arr[,,l]%*%model$fhatk[,l]))
          dzdfhatl <- model$y[,l]-exp(fhat)
          #first
          output[i,j] <- output[i,j] + as.numeric(dzdfhatl%*%dfhatl)
          #second
          output[i,j] <- output[i,j] - 0.5*(t(dfhatl)%*%Kxinv%*%fhat-t(fhat)%*%Kxinv%*%dKxdijmat%*%Kxinv%*%fhat+t(fhat)%*%Kxinv%*%dfhatl)
          
          #third
          #first3 <- 0.5*fhat%*%Kxinv%*%dKxdijmat%*%Kxinv%*%fhat
          prod <- dKxdijmat%*%model$Skinv_arr[,,l]
          
          second3 <- 0.5*sum(diag(inv%*%prod))
          output[i,j] <- output[i,j]-second3
          #output[i,j] <- output[i,j] + 
        }
      }
    }
    return(as.vector(-output+model$Ktinv%*%xmat))
  }
  output <- scg(objfunc,gradObj,as.vector(model$x),1)
  return(output)
}
