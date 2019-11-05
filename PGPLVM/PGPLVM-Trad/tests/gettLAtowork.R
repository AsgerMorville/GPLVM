#Project: Get tLA to work

sourceDirectory("Modules/",modifiedOnly = F)

model1 <- generateData(30,50,1)

plot(model1$x)

#Fit it with the tla
Xny <- eigenInit(model1$y,1)


model <- list()
model$x <- Xny
model$q <- dim(model$x)[2]
model$y <- model1$y
model$p <- dim(model1$y)[2]
model$t <- dim(model$x)[1]

model$maxiter <- 100
model$gradopt <- "tLa"
Kt <- fillup(as.matrix(1:model$t),ARDKernel,model1$thetax)
model$Ktinv <- solve(Kt)


vecLogLik(as.vector(Xny))
vecGrad(as.vector(Xny))
end <- scg(vecLogLik,vecGrad,as.vector(Xny),1)
plot(end)
plot(model1$x)
##
model$gradopt <- "aLa"
end2 <- scg(vecLogLik,vecGrad,as.vector(Xny),1)
plot(model1$x)
plot(end2)
