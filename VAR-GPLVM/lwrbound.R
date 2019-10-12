library(R.utils)
sourceDirectory("timedyn/Modules/",modifiedOnly = F)

thetaf=list(sigmaf=1,lengthscales=1)
theta = thetaf
n <- t <- 30
q <- 1
m<-30
p <- 20
#generate data object
data <- generateData(t,p,q,parameters=theta)
#data$x <- scale(data$x,center=T,scale=F)

plot(1:30,data$x)

#Test of implementation
indices <- seq(from=1,to=n,by=1)

modeltest <- list()
modeltest$tvec <- as.matrix(1:30)
modeltest$y <- data$y
modeltest$xu <- as.matrix(data$x[indices])

#lambda
modeltest$thetaf <- list(sigmaf=1,lengthscales=1)
modeltest$thetax <- list(l=1,r=1)
modeltest$beta <- 100


#mustart <- as.matrix(data$x)
mustart <- (data$kx)%*%(1000*eigenInit(data$y,1))

#calibrateLambda(n,q,as.matrix(1:n))
#Sjcalculator(matrix(runif(n*q,min=0.4,max=0.6)),as.matrix(1:n),n,q)
obj <- calibrateLambda(n,q,modeltest$tvec)
#Lstart <- obj$covars
#modeltest$Lambda <- obj$covars-0.5

#Check
Sjcalculator(obj$covars,as.matrix(1:30),n,q)
Lstart <- matrix(obj$covars,nrow=t,ncol=q)
#xustart <- as.matrix(mustart[indices])
#xustart <- as.matrix(runif(5))
plot(1:30,data$x)
plot(1:30,mustart)

#startguess <- listToVec(list(mubar=mustart,L=Lstart,xu=xustart))
startguess <- listToVec(list(mubar=mustart,L=Lstart))
#startguess <- listToVec(list(mubar=data$x,L=Lstart))
#startguess<- mustart
list <- vecToList(startguess,t,q,p,m)


#plot(1:100,-mustart)

#end <- scg(vecLogLik,vecGradLogLik,startguess,1)
end <- constrained_scg(vecLogLik,vecGradLogLik,startguess,1,c(31:60))
vecLogLik(end)
vecLogLik(startguess)
plot(1:t,data$x,col="red")
plot(1:t,-end[1:t])
plot(1:t,-mustart)
plot(1:t,(data$kx)%*%end[1:t])

plot(1:t,-end[1:t])
plot(1:t,end[(t+1):(2*t)])
plot(1:t,startguess[(t+1):(2*t)])
plot(1:60,vecGradLogLik(startguess))
plot(1:60,vecGradLogLik(end))
modeltest$mubar <-as.matrix(end[1:100])
modeltest$Lambda <- as.matrix(end[101:200])
#modeltest$xu <- as.matrix(end[201:220])
kx <-modCalculator(modeltest)$Kx
kx-data$kx
plot(1:100,kx%*%end[1:100])

Sjcalculator(matrix(1,ncol=q,nrow=n),as.matrix(1:30),n,q)



plot(1:100,data$x)


plot(1:100,data$x)

modeltest$xu


