#Oil data
data <- read.table("DataTrn.txt")
lbl <- read.table("DataTrnLbls.txt")

Ysc <- scale(data,scale=F)
YYt <- Ysc%*%t(Ysc)

alpha <- 1
beta <- 1

eig <- eigen(YYt)
D <- 12
l1 <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
l2 <-1/sqrt(eig$values[2]/(alpha*D)-1/(alpha*beta))



X1 <- eig$vectors[,1]*l1
X2 <- eig$vectors[,2]*l2

plot(X1,X2)
plot(Xny,y[,2])
lblvec <- numeric(length(X1))
for (i in 1:length(X1)){
  if (lbl[i,1]==1){
    lblvec[i] <- "A"
  }
  if (lbl[i,2]==1){
    lblvec[i] <- "B"
  }
  if (lbl[i,3]==1){
    lblvec[i] <- "C"
  }
}
lblvec
library(ggplot2)

df <- data.frame(cbind(X1,X2))
df$xx <- as.factor(lblvec)
#df$X1 <- as.numeric(df$X1)
#df$X2 <- as.numeric(df$X2)
ggplot(data=df)+geom_point(data=df,aes(x=X1,y=X2,color=xx))
## Sub sample

subdata <- data[1:100,]
sublbl <- lbl[1:100,]

Ysc <- scale(subdata,scale=F)
YYt <- Ysc%*%t(Ysc)

alpha <- 1
beta <- 1

eig <- eigen(YYt)
D <- 12
l1 <- 1/sqrt(eig$values[1]/(alpha*D)-1/(alpha*beta))
l2 <-1/sqrt(eig$values[2]/(alpha*D)-1/(alpha*beta))



X1 <- eig$vectors[,1]*l1
X2 <- eig$vectors[,2]*l2

plot(X1,X2)
plot(Xny,y[,2])
lblvec <- numeric(length(X1))
for (i in 1:length(X1)){
  if (sublbl[i,1]==1){
    lblvec[i] <- "A"
  }
  if (sublbl[i,2]==1){
    lblvec[i] <- "B"
  }
  if (sublbl[i,3]==1){
    lblvec[i] <- "C"
  }
}
lblvec
library(ggplot2)

df <- data.frame(cbind(X1,X2))
df$xx <- as.factor(lblvec)
#df$X1 <- as.numeric(df$X1)
#df$X2 <- as.numeric(df$X2)
ggplot(data=df)+geom_point(data=df,aes(x=X1,y=X2,color=xx))
