#Optim test

n <- 100

x <- seq(from=-5,to=5,length.out = n)



Kernel <- fillup(x,x,kern,3)

Y <- mvrnorm(n=1,mu=rep(0,n),Sigma = Kernel)
Y2 <- Y+rnorm(n,sd=0.1)

plot(x,Y2)
lines(x,Y)

y <- Y2
optim(c(1),simpelloglik,simpelGrad,method ="BFGS" )

,lower=c(0.01),method="L-BFGS-B"
simpelloglik(3)

Kernel <- fillup(x,x,kern,3)
invx <- solve(Kernel+diag(nrow=n)*0.1^2)
(-0.5*t(y)%*%invx%*%y-0.5*log(abs(det(Kernel+diag(nrow=n)*0.1^2)))

0.5*log(abs(det(Kernel+diag(nrow=n)*0.1^2)))
det(Kernel+diag(nrow=n)*0.1)
simpelloglik(3)
gridl <- seq(from=0.01,to=50,length.out=100)
count <- 1
tester <- numeric(100)
for (i in gridl){
  tester[count] <- simpelGrad(i)
  count<- count+1
}
plot(tester)
simpelloglik(c(1,2))
