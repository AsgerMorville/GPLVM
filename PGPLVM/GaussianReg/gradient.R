

gradDescent <- function(grad,init,gamma){
  a <- init
  for (i in 1:100){
    nab <- grad(a)
    a <- a-gamma*nab
  }
  return(a)
}

gradd <- function(x){
  return(4*x-2)
}

gradDescent(gradd,0,0.1)

