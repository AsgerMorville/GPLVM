library(R.utils)
sourceDirectory("timedyn/Modules/",modifiedOnly = F)

thetaf=list(sigmaf=1,lengthscales=1)
theta = thetaf
n <- 30
t <- 30
q <- 1
m<-30
p <- 20
#generate data object
data <- generateData(t,p,q,parameters=theta)
#data$x <- scale(data$x,center=T,scale=F)

plot(1:30,data$x)

#Test of implementation
indices <- floor(seq(from=1,to=n,length.out = m))

modeltest <- list()
modeltest$tvec <- as.matrix(1:30)
modeltest$y <- data$y
#modeltest$xu <- as.matrix(data$x[indices])

#lambda
modeltest$thetaf <- list(sigmaf=1)
#modeltest$thetaf <- thetaf
modeltest$thetax <- list(l=1,r=1)
modeltest$beta <- 50


#mustart <- as.matrix(data$x)
mustart <- solve(data$kx)%*%(100*eigenInit(data$y,1))

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
plot(1:30,-mustart)

xustart <- (100*eigenInit(data$y,1))[indices]
#startguess <- listToVec(list(mubar=mustart,L=Lstart,xu=xustart))
maxlat <- apply(mustart,2,max)
minlat <- apply(mustart,2,min)
lengthscalestart <- 1/(maxlat-minlat)^2
startguess <- listToVec(list(mubar=mustart,L=Lstart,xu=as.matrix(data$x[indices]),lengthscales=0.17))



#end <- scg(vecLogLik,vecGradLogLik,startguess,1)
end <- constrained_scg(vecLogLik,vecGradLogLik,startguess,1,c(31:60))




library(numDeriv)
num <- grad(vecLogLik,startguess)
reel <- vecGradLogLik(startguess)
plot(num-reel)

#startguess <- listToVec(list(mubar=data$x,L=Lstart))
#startguess<- mustart
list <- vecToList(startguess,t,q,p,m)


vecLogLik(startguess)
vecGradLogLik(startguess)
#plot(1:100,-mustart)

#Numerical optimizer
numgrad <- function(guess){
  return(grad(vecLogLik,guess))
}


plot(end-startguess)
vecLogLik(end)
vecLogLik(startguess)
plot(end-startguess)
num <- grad(vecLogLik,startguess)
reel <- vecGradLogLik(startguess)
plot(num)
plot(reel)
plot(reel-num)
plot(vecGradLogLik(startguess))
plot(jacobian(vecGradLogLik,startguess))
plot(1:t,data$x,col="red")
plot(1:t,-(data$kx)%*%end[1:t])
plot(1:t,(data$kx)%*%mustart)
plot(1:t,(data$kx)%*%end[(t+1):t])
plot(mustart,end[1:t])
#Test of gradients wrt xu
#first 
modeltest$mubar <-mustart
modeltest$Lambda <- Lstart
plot(startguess-end)
plot(end)

testfunc <- function(x){
  modeltest$xu <- as.matrix(x)
  mod1 <- modCalculator(modeltest)
  return(dPsi2dxu(1,1,mod1)[1,1])
}

testfunc()
plot(1:t,end[1:t])
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


