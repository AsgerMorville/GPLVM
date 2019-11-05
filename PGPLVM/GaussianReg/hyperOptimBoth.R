#HyperOptim, sigma and length scale

n <- 100

x <- seq(from=-5,to=5,length.out = n)



Kernel <- fillup(x,x,kern,1)
 
Y <- mvrnorm(n=1,mu=rep(0,n),Sigma = Kernel)
Y2 <- Y+rnorm(n,sd=0.1)

plot(x,Y2)
lines(x,Y)

y <- Y2
optim(c(0.5,0.5),loglik,grad5_9,method="BFGS")


points<-Cholpred(c(3,0.1),seq(from=-5,to=5,length.out = n))
points <- Cholpred(c(0.3099558,-0.1314527),seq(from=-5,to=5,length.out = n))
points <- Cholpred(c(1,-0.1314527),seq(from=-5,to=5,length.out = n))
plot(seq(from=-5,to=5,length.out=n),points)

 
#Plot
matrixx <- matrix(NA,nrow=100,ncol=100)
gridl <- seq(from=0.1,to=30,length.out=100)
grids <- seq(from=0.003,to=0.1,length.out=n2)

for (i in 1:100){
  for (j in 1:100){
    matrixx[i,j] <- loglik(c(gridl[i],grids[j]))
  }
  print(i)
}

image(matrixx,asp=1)
matrixx[50,3]
loglik(c(2,0.003))
matrixx[20,10]
max(matrixx)
gridl[20]
grids[10]
n2 <- 10000
grids
tester <- numeric(n2)
for(i in 1:n2){
  tester[i] <- norm(grad5_9(c(5.98,grids[i])),type="2")
  print(i)
}
plot(tester)
tester
min(tester)
