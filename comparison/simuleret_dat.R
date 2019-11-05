#First create the object that we want to keep throughout the session. Call this for MODEL0.
#Then use rm(list=setdiff(ls(), "MODEL0")) to keep it when clearing the sessions for new GPLVM model

library(R.utils)
library(MASS)
setwd('C:/Users/Asger/Google Drive/Speciale/')

#First, genereate data
sourceDirectory("PGPLVM/dLA/Modules/",modifiedOnly = F)

#generate the data
MODEL0 <- generateData(t=30,p=100,q=1)
MODEL0$Xny <-1000*eigenInit(MODEL0$y,1)

#dLA
dataobj <- list()
dataobj$y <- MODEL0$y
dataobj$x0 <- MODEL0$Xny
dataobj$maxiter <- 100
dataobj$thetaf <- MODEL0$thetaf
dataobj$thetax <- MODEL0$thetax
bbj <- dLA(dataobj)

MODEL0$dLAobj <- bbj

plot(MODEL0$x)
lines(-MODEL0$dLAobj)
#plot(MODEL0$dLAobj)


rm(list=setdiff(ls(), "MODEL0"))
#tLA
sourceDirectory("PGPLVM/PGPLVM-Trad/Modules/",modifiedOnly = F)
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

end <- scg(vecLogLik,vecGrad,as.vector(model$x),1)

MODEL0$tLA <- end

#aLa

model$gradopt <- "aLa"
MODEL0$aLA <- scg(vecLogLik,vecGrad,as.vector(model$x),1)

plot(MODEL0$x)
lines(MODEL0$dLAobj,col="red")
lines(MODEL0$tLA,col="blue")
lines(MODEL0$aLA,col="green")

rm(list=setdiff(ls(), "MODEL0"))

#Try with GPLVM
sourceDirectory("GaussianReg/GPLVM-trad/Modules/",modifiedOnly=FALSE)
tvec <- seq(from=0,to=MODEL0$t,length.out = MODEL0$t)
Kt <- fillup(tvec,tvec,matern,1)

modeltest <- list()
#modeltest$x <- x
modeltest$y <- MODEL0$y
modeltest$hyper$gamma <- 1
modeltest$hyper$alpha <- 1
modeltest$YYt <- MODEL0$y%*%t(MODEL0$y)
modeltest$Ktinv <- solve(Kt)
#modeltest$hyper$beta <- 100
modeltest$q <- dim(MODEL0$Xny)[2]

vecLogLik(c(1,1,100,as.vector(MODEL0$Xny)))
vecGrad()
GPLVM <- scg(vecLogLik,vecGrad,c(1,1,100,as.vector(MODEL0$Xny)),3)
MODEL0$GPLVM <- GPLVM

#How about Var-GPLVM

gplmpoi <- GPLVM[-c(1:3)]
plot(MODEL0$x,ylim=c(-3,3))
lines(-MODEL0$dLAobj,col="red")
lines(-MODEL0$tLA,col="blue")
lines(-MODEL0$aLA,col="green")
lines(-gplmpoi/1.5)

lines(MODEL0$GPLVM/1.,col="black")
plot(MODEL0$x)
lines(GPLVM[-c(1:3)])

library(numDeriv)
num <- grad(vecLogLik,c(3,1.67,30,as.vector(MODEL0$Xny)))
reel <- vecGrad(c(3,1.67,30,as.vector(MODEL0$Xny)))
plot(num-reel)
num
reel
