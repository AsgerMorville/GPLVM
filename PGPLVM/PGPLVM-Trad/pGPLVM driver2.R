#Approximated Laplace 
set.seed(3)
library(MASS)

source("Modules/funks.R")
model <- list()
model$x <- MODEL0$Xny
model$q <- dim(model$x)[2]
model$p <- dim(MODEL0$y)[2]
model$t <- dim(model$x)[1]
model$y <- MODEL0$y
model$maxiter <- 100
model$gradopt <- "tLa"
Kt <- fillup(as.matrix(1:model$t),ARDKernel,MODEL0$thetax)
model$Ktinv <- solve(Kt)

vecLogLik(as.vector(model$x))
vecGrad(as.vector(model$x))

end <- scg(vecLogLik,vecGrad,as.vector(model$x),1)


plot(MODEL0$x)
plot(-end)
plot(-MODEL0$Xny)

testlist$x<- MODEL0$Xny
testlist$y <- MODEL0$y
testlist$thetaf <- MODEL0$thetaf
testlist$thetax <- MODEL0$thetax
tLA(testlist)


tLA
vecGrad(as.vector(MODEL0$x))

test<-scg(testfunction,testgrad,Xny,1)
test2 <- scg(vecLogLik,vecGrad,Xny,1)
plot(1:100,-Xny)
plot(1:100,x)
plot(1:100,-test2)
plot(tvec,x)
plot(tvec,test)
plot(x,f[2,])
plot(-realx2,f[2,])
plot(-realx2,f[2,],xlim=c(-3,3))
plot(x,f[2,],xlim = c(-3,3))