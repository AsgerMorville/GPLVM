#dLA driver
library(R.utils)
library(MASS)
sourceDirectory("Modules/",modifiedOnly = F)

#generate the data
model1 <- generateData(t=30,p=60,q=1)

#plot the data
for (h in 1:model1$p){
  plot(model1$y[,h])
  lines(exp(model1$f[,h]))
  Sys.sleep(1)
}

#find the eigenvalue init
Xny<-1000*eigenInit(model1$y,1)

#Compare with the true x

plot(Xny[,1])
plot(model1$x[,1])

#Use the dLA to find a better estimate of the latent x-matrix.
#the dLA takes a list containing observations and start guess of latent matrix

dataobj <- list()
dataobj$y <- model1$y
dataobj$x0 <- Xny
dataobj$maxiter <- 100
dataobj$thetaf <- model1$thetaf
dataobj$thetax <- model1$thetax
obj <- dLA(dataobj)

plot(model1$x[,1])
plot(-obj)

plot(model1$x[1:50,2])
plot(-obj[101:150])
