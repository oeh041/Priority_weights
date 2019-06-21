## Priority weights
setwd("yourfolder")
if(F){
  rm(list=ls())
  gc()
  cat("\014")
}
## w = a + (age-b)*c*exp(-(age-b)*d)
## anchor: age at which weight=a
## a: weight at age=anchor (up-down shift)
## b: right-left shift
## c: scale parameter
## d: rate at which the weight approaches a
# weight=function(age,a=1,b=0,c=0.165,d=0.06,anchor=80){
#   w=c()
#   w=a+(age-b)*c*exp(-(age-b)*d)
#   w=w-w[age==anchor]+a
#   w[age>anchor]=a
#   w
# }
age <- 0:100
## Anchors at a for age=anchor, and then sets weight=a
weight <- function(age,anchor=80,a=1,b=0,c=0.165,d=0.06){
  if(sum(age<0)>0) warning("Negative ages")
  age <- age[age>=0]
  w <- a+c*((age-b)*exp(-(age-b)*d)-(anchor-b)*exp(-(anchor-b)*d))
  w[age>anchor] <- a
  w
}
## weight2 sums to maxage. weight2==1 for all ages corresponds to uniform distribution.
weight2 <- function(age,maxage=100,b=0,c=0.165,d=0.06){
  a <- (maxage*d-(1-exp(maxage*d))*(b*d-1))*c*exp(-(maxage-b)*d)/(maxage*d^2)+1
  age <- age[age>=0];if(sum(age<0)>0) warning("Negative ages")
  w <- a+c*(age-b)*exp(-(age-b)*d)
  w
}
## Simplest form
weight3 <- function(a,b,c,d,age){
  a+(age-b)*c*exp(-d*(age-b))
}
## Anchors at a for age=anchor
weight4 <- function(age,anchor=80,a=1,b=0,c=0.165,d=0.04){
  if(sum(age<0)>0) warning("Negative ages")
  age <- age[age>=0]
  w <- a+c*((age-b)*exp(-(age-b)*d)-(anchor-b)*exp(-(anchor-b)*d))
  w
}

## Least squares
## Minimizes SUM{(PW_i-(a-(t_i-d)*c*exp(-b(t_i-d))))^2}
##          =SUM{(PW_i-PW.modell_i)^2}
## Finds minimum for a, b, c og d
## Use simplest form (weight3)
## The data
spm8a <- as.matrix(read.csv2("Spm8a.csv"))
head(spm8a);class(spm8a);class(spm8a[,1])
spm8a <- na.omit(spm8a);dim(spm8a)
sum(spm8a[,"konsistens"])
spm8a <- spm8a[spm8a[,"konsistens"]==1,];dim(spm8a);head(spm8a)
inds <- nrow(spm8a)
X <- matrix(NA,nrow=inds,ncol=5);colnames(X)=c(10,25,40,55,70);head(X)
X[,1:4] <- 10/spm8a[,1:4];head(X)
X[,5] <- 1;head(X)

## Estimating all parameters without restrictions
## Includes data for reference category PW(70)=1
pw.ls <- function(x){
  ts <- as.numeric(rep(colnames(X),each=inds));ts
  pws <- as.numeric(X);pws
  sum((pws-weight3(a=x[1],b=x[2],c=x[3],d=x[4],ts))^2)
}
means <- optim(par=rep(0,4),fn=pw.ls)$par;means
## Estimating all parameters without restrictions
## Does NOT include data for reference category PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum((pws-weight3(a=x[1],b=x[2],c=x[3],d=x[4],ts))^2)
}
means2 <- optim(par=rep(0,4),fn=pw.ls)$par;means2
## Fixing the max at t=10
## Does NOT include data for reference category PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum((pws-weight3(a=x[1],b=10-1/x[3],c=x[2],d=x[3],ts))^2)
}
means3 <- optim(par=rep(1,3),fn=pw.ls)$par;means3
## Fixing the max at t=10
## Include data for reference category PW(70)=1
pw.ls <- function(x){
  ts <- as.numeric(rep(colnames(X),each=inds));ts
  pws <- as.numeric(X);pws
  sum((pws-weight3(a=x[1],b=10-1/x[3],c=x[2],d=x[3],ts))^2)
}
means4 <- optim(par=rep(1,3),fn=pw.ls)$par;means4
## Does NOT include data for reference category PW(70)=1
## Anchoring at PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum((pws-weight4(anchor=70,a=1,b=x[1],c=x[2],d=x[3],ts))^2)
}
means5 <- optim(par=rep(1,3),fn=pw.ls)$par;means5
## Fixing the max at t=10
## Include data for reference category PW(70)=1
## Anchoring at PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum((pws-weight4(anchor=70,a=1,b=10-1/x[2],c=x[1],d=x[2],ts))^2)
}
means6 <- optim(par=rep(1,2),fn=pw.ls)$par;means6
## Coefficients
## gamma
b <- 10-1/means6[2];b
## C
c <- means6[1];c
## beta
d <- means6[2];d
## alpha
a <- 1-(70-b)*c*exp(-d*(70-b));a
weight3(a = a, b = b, c = c, d = d, age = 10:70)

Xmeans <- X[,c("10","25","40","55")];Xmeans
ts.means <- as.numeric(rep(colnames(Xmeans),each=inds));ts.means
pws.means <- as.numeric(Xmeans);pws.means
res.means <- pws.means-weight4(anchor=70,a=1,b=10-1/means6[2],c=means6[1],d=means6[2],age=ts.means);round(res.means,2)

## MEDIAN
## Estimating all parameters without restrictions
## Includes data for reference category PW(70)=1
pw.ls <- function(x){
  ts <- as.numeric(rep(colnames(X),each=inds));ts
  pws <- as.numeric(X);pws
  sum(abs(pws-weight3(a=x[1],b=x[2],c=x[3],d=x[4],ts)))
}
medians <- optim(par=rep(0,4),fn=pw.ls)$par;medians
## Estimating all parameters without restrictions
## Does NOT include data for reference category PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum(abs(pws-weight3(a=x[1],b=x[2],c=x[3],d=x[4],ts)))
}
medians2 <- optim(par=rep(0,4),fn=pw.ls)$par;medians2
## Fixing the max at t=10
## Does NOT include data for reference category PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum(abs(pws-weight3(a=x[1],b=10-1/x[3],c=x[2],d=x[3],ts)))
}
medians3 <- optim(par=rep(1,3),fn=pw.ls)$par;medians3
## Fixing the max at t=10
## Include data for reference category PW(70)=1
pw.ls <- function(x){
  ts <- as.numeric(rep(colnames(X),each=inds));ts
  pws <- as.numeric(X);pws
  sum(abs(pws-weight3(a=x[1],b=10-1/x[3],c=x[2],d=x[3],ts)))
}
medians4 <- optim(par=rep(1,3),fn=pw.ls)$par;medians4
## Does NOT include data for reference category PW(70)=1
## Anchoring at PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum(abs(pws-weight4(anchor=70,a=1,b=x[1],c=x[2],d=x[3],ts)))
}
medians5 <- optim(par=rep(1,3),fn=pw.ls)$par;medians5
## Fixing the max at t=10
## Include data for reference category PW(70)=1
## Anchoring at PW(70)=1
pw.ls <- function(x){
  X2 <- X[,c("10","25","40","55")];X2
  ts <- as.numeric(rep(colnames(X2),each=inds));ts
  pws <- as.numeric(X2);pws
  sum(abs(pws-weight4(anchor=70,a=1,b=10-1/x[2],c=x[1],d=x[2],ts)))
}
medians6 <- optim(par=rep(1,2),fn=pw.ls)$par;medians6
## Coefficients
## gamma
b <- 10-1/medians6[2];b
## C
c <- medians6[1];c
## beta
d <- medians6[2];d
## alpha
a <- 1-(70-b)*c*exp(-d*(70-b));a
weight3(a = a, b = b, c = c, d = d, age = 10:70)

Xmedians <- X[,c("10","25","40","55")];Xmedians
ts.medians <- as.numeric(rep(colnames(Xmedians),each=inds));ts.medians
pws.medians <- as.numeric(Xmedians);pws.medians
res.medians <- pws.medians-weight4(anchor=70,a=1,b=10-1/medians6[2],c=medians6[1],d=medians6[2],age=ts.medians);round(res.medians,2)

## MEAN (red/black +) and MEDIAN (blue/magenta x)
tiff(filename = "FigWeight.tiff",units = "in",width = 10, height = 10,
     compression = "lzw",res=300)
# windows(10,10)
par(mfrow=c(1,1))
plot(NA,ylim=range(X),xlim=range(as.integer(colnames(X))),xlab="QALE",ylab="Weight",cex.lab=1.5,cex.axis=1.5)
for(i in 1:inds){for(j in colnames(X)[-5]){
  points(as.integer(j),X[i,j],cex=2*sqrt(sum(X[,j]==X[i,j])),col="grey",pch=16)}
}
abline(h=c(seq(0,5,1),seq(10,20,5)),lty=3)
points(colnames(X)[-5],apply(X,2,mean)[-5],pch=3,cex=2,lwd=4)
points(colnames(X)[-5],apply(X,2,median)[-5],pch=4,cex=2,lwd=4,col="blue")
points(70,1,pch=16,col="black")
text(70,1,"ref",pos=3,cex=1.5)
lines(10:70,weight4(anchor=70,a=1,b=10-1/means6[2],c=means6[1],d=means6[2],age=10:70),lwd=2,col="red")
lines(10:70,weight4(anchor=70,a=1,b=10-1/medians6[2],c=medians6[1],d=medians6[2],age=10:70),lwd=2,col="magenta")
legend(x=10,y=19.8,pch=c(NA,NA,3,4,16),col=c("red","magenta","black","blue","grey"),
       legend=c("Estimated line, mean","Estimated line, median","Observed means","Observed medians","All observations, scaled by frequency"),
       lty=c(1,1,NA,NA,NA),cex=1.5,bg="white",box.col="white",lwd=3)
dev.off()
shell("open FigWeight.tiff")
