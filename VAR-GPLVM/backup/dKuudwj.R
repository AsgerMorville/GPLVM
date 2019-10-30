dKuudwj <- function(j,mod){
  Kuu <- mod$Kuu
  xu <- mod$xu[,j]
  factor <- outer(xu,xu,'-')^2
  return(Kuu*factor/(-2))
}
