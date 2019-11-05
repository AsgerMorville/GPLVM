optimGPLVM <- function(Y,dim,startpar,xstart,maxiter){
  n <- dim(Y)[1]
  #Xny <- init(Y,dim)
  par <- startpar
  xx <- xstart
  for (i in 1:maxiter){
    #Parameter update
    parfunc <- function(par){
      return(loglikchol(c(par,xx)))
    }
    pargrad <- function(par){
      return(gradientChol2.0(c(par,xx))[1:3])
    }
    par <- scg(parfunc,pargrad, par,0.1)
    
    #x update 
    for (i in 1:n){
      xfunc <- function(xxx){
        if (i<n && i >1){
          return(loglikchol(c(par,xx[1:(i-1)],xxx,xx[(i+1):n])))
        } 
        if (i == 1){
          return(loglikchol(c(par,xxx,xx[(i+1):n])))
        }
        if (i == n){
          return(loglikchol(c(par,xx[1:(i-1)],xxx)))
        }
      }
      xgrad <- function(xxx){
        if (i<n && i >1){
          return(gradientChol2.0(c(par,xx[1:(i-1)],xxx,xx[(i+1):n]))[3+i])
        } 
        if (i == 1){
          return(gradientChol2.0(c(par,xxx,xx[(i+1):n]))[3+i])
        }
        if (i == n){
          return(gradientChol2.0(c(par,xx[1:(i-1)],xxx))[3+i])
        }
      }
      xx[i] <- scg(xfunc,xgrad,xx[i],0.1)
      
      
    }
    if (norm(gradientChol2.0(c(par,xx)),type="2")<1){break}
  }
  return(c(par,xx))
}
