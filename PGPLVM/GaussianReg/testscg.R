testfunc <- function(x){
  return((1-x[1])^2+100*(x[2]-x[1]^2)^2)
}
testgrad <- function(x){
  return(c(-400*x[1]*(x[2]-x[1]^2)-2*(1-x[1]),200*(x[2]-x[1]^2)))
}

source("Scgfrommatlab.R")
testfunc(c(1,1.2))
testgrad(c(1,1.2))
scg(testfunc,testgrad,c(-1.2,1),0.0001)

par <-optim(c(-1.2,1),testfunc,method="CG")
par
paste(hej$rCode)


hej <- mat2r(matIn, verbose = 0)$rCode
hej <- mat2r(matIn, verbose = 0)
hej$rCode

testgrad(par$par)
?optim
scg <- function(init,fn,gn,sig,lam1,maxiter){
  w <- init
  lam1bar <- 0
  p <- -gn(w1)
  r1 <- p
  k <- 1
  success <- TRUE
  sigmas <- numeric(maxiter)
  
  for(k in 1:maxiter){
    ## 2
    if (success == TRUE){
      sigmas[k] <- sig/norm(p,type="2")
      s <- (gn(w+sigmas[k]*p)-gn(w))/sigmas[k]
      delta <- t(p)%*%s
    }
    ## 3 scale sk
    
    
    
  }
}