#Chol inverse
chol_inverse <- function(A){
  L <- t(chol(A))
  Linv <- solve(L)
  return(list(Kinv=t(Linv)%*%Linv,L=L))
}
