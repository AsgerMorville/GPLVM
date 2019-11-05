hessianNewton <- function(f,Kx_inv){
  return(diag(-exp(f))-Kx_inv)
}