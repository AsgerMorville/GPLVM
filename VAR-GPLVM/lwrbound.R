library(R.utils)
sourceDirectory("Modules/",modifiedOnly = F)

thetaf=list(sigmaf=1,lengthscales=1)
theta = thetaf
n <- t <- 100
q <- 1
m<-50
p <- 5
#generate data object
data <- generateData(t,p,q,parameters=theta)


plot(1:100,data$x)

#Test of implementation
indices <- seq(from=1,to=100,by=2)

modeltest <- list()
modeltest$tvec <- as.matrix(1:100)
modeltest$y <- data$y
#modeltest$xu <- as.matrix(data$x[indices])

#lambda
modeltest$thetaf <- list(sigmaf=1,lengthscales=1)
modeltest$thetax <- list(l=20,r=1)
modeltest$beta <- 10

#mustart <- as.matrix(data$x)
mustart <- eigenInit(data$y,1)
Lstart <- matrix(rnorm(t*q,mean=0.5,sd=0.1),nrow=t,ncol=q)
xustart <- as.matrix(mustart[indices])


startguess <- listToVec(list(mu=mustart,L=Lstart,xu=xustart))
plot(1:100,mustart)

end <- scg(vecLogLik,vecGradLogLik,startguess,1)
plot(1:100,mustart)
plot(1:100,data$x)

plot(1:100,end[1:100])


modeltest$xu








n<- 5
q <- 3 
p <- 4
m <-2
modd <- list()
modd$y <- matrix(1:(n*p),nrow=n,ncol=p,byrow=F)
modd$mubar <- matrix(1:(n*q),nrow=n,byrow=F)
modd$Lambda <- matrix(runif(n*q),nrow=n,ncol=q)
modd$tvec <- as.matrix(1:n)
modd$xu <- matrix(1:(m*q),nrow=m,ncol=q,byrow=F)
modd$thetaf <- list(w=1:q,sigmaf=1)
modd$thetax <- list(r=1,l=1)
modd$beta <- 100
#modd$hyper <- list(w=1:q,sigmaf=1,sigma=1,par=list(r=1,l=1))


#Main script: we are given data Y in R t times p. We decide latent dim, q. 

#ModInit: initial guess
ModInit <- list()


mod1 <- modCalculator(modd)

logLik2(modd)
dFhatdmu(mod1)



S1 <- matrix(NA,ncol=q,nrow=n)
for (i in 1:q){
  S1[,i]<- diag(as.matrix(S[,,q]))
}
psistatsobj <- psistats(mu,S,xu,hyper,n,m)
psistatsobj

dPsi1dmu(2,3,psistatsobj$psi1,hyper$w,mu,S1,xu)
dPsi2dmu(2,3,psistatsobj$psi2,hyper$w,mu,S1,xu)


logLik(tvec,y,mu,S,xu,hyper)
testmat <- matrix(c(2,0,0,1),ncol=2)
