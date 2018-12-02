library(readr)
stars <- read_csv("/SDSS_stars.csv")
pairs(stars)
NS=stars=Neutron.Stars

plot(NS[,1:2],col=NS$ng)


NS1=subset(NS,NS$ng==1)
plot(NS1)

library(np)
M0=lm(NS1$lb~NS1$lp)
summary(M0)
PM0=predict(M0,interval = 'confidence')

RT=c("lc","ll")
BWM=c("cv.aic","cv.ls")
BWT=c("fixed","generalized_nn","adaptive_nn")
KT=c("gaussian","epanechnikov","uniform")

n1=length(RT)
n2=length(BWM)
n3=length(BWT)
n4=length(KT)


T=matrix(nrow=n1*n2*n3*n4,ncol=8)
k=1
for(i1 in 1:n1)
{
  for(i2 in 1:n2)
  {
    for(i3 in 1:n3)
    {
      for(i4 in 1:n4)
      {
        bw=npregbw(NS1$lb~NS1$lp,regtype=RT[i1],bwmethod=BWM[i2],bwtype=BWT[i3],ckertype=KT[i4])
        R=npreg(bw)
        T[k,1:4]=c(RT[i1],BWM[i2],BWT[i3],KT[i4])
        T[k,5:8]=c(R$R2,R$MSE,R$MAE,R$MAPE)
        k=k+1
      }
    }
  }
}



OT=T[order(T[,5],decreasing = TRUE),]
colnames(OT)=c("RegType","BWMethod","BWType","KernType","R^2","MSE","MAE","MAPE")
View(OT)


bw=npregbw(NS1$lb~NS1$lp,regtype=OT[1,1],bwmethod=OT[1,2],bwtype=OT[1,3],ckertype=OT[1,4])

R=npreg(bw)
plot(NS1[,1:2],ylim=c(7,16),col="gray75")
abline(M0,col=2)
lines(unlist(NS1[,1]),PM0[,2],lty=2)
lines(unlist(NS1[,1]),PM0[,3],lty=2)
par(new=TRUE)
plot(R, plot.errors.method="asymptotic",col=2,ylim=c(7,16),xlab="",ylab="")
par(new=TRUE)
plot(R, plot.errors.method="bootstrap",col=2,ylim=c(7,16),xlab="",ylab="")


bw=npregbw(NS1$lb~NS1$lp,regtype=OT[1,1],bwmethod=OT[1,2],bwtype=OT[1,3],ckertype=OT[2,4])



###############


fhat <- npcdens(NS1$lb~NS1$lp, tol = 0.1, ftol = 0.1)
fhatunc<-npudens(~NS1$lb)
plot(fhat)
plot(fhatunc)
names(fhat)

CD=cbind(unlist(fhat$yeval),fhat$condens)
CD=CD[order(CD[,1]),]
plot(CD,type="l",xlab="y",ylab="f(y)")

UD=cbind(unlist(fhatunc$eval),fhatunc$dens)
UD=UD[order(UD[,1]),]
lines(UD,type="l",col=2)
legend("topright",legend=c("Conditional Density","Unconditional Density"),col=1:2,lty=1)



###########################

fhat <- npcdist(NS1$lb~NS1$lp)
names(fhat)
fhatunc<-npudist(~NS1$lb)
plot(fhat)
plot(fhatunc)
names(fhat)

CD=cbind(unlist(fhat$yeval),fhat$condist)
CD=CD[order(CD[,1]),]
plot(CD,type="l",xlab="y",ylab="f(y)")

UD=cbind(unlist(fhatunc$eval),fhatunc$dist)
UD=UD[order(UD[,1]),]
lines(UD,type="l",col=2)
legend("bottomright",legend=c("Conditional Density","Unconditional Density"),col=1:2,lty=1)

############################


bw <- npcdistbw(NS1$lb~NS1$lp,tol=0.01,ftol=0.01)
model.q0.05 <- npqreg(bws = bw, tau = 0.05)
model.q0.25 <- npqreg(bws = bw, tau = 0.25)
model.q0.50 <- npqreg(bws = bw, tau = 0.50)
model.q0.75 <- npqreg(bws = bw, tau = 0.75)
model.q0.95 <- npqreg(bws = bw, tau = 0.95)
NSD=cbind(NS1$lp,NS1$lb,model.q0.05$quantile,model.q0.25$quantile,model.q0.50$quantile,model.q0.75$quantile,model.q0.95$quantile)
ONSD=NSD[order(NSD[,1]),]
plot(ONSD[,1],ONSD[,2],main = "",xlab = "Period", ylab = "Magnetic Field",col="gray75")
lines(ONSD[,1],ONSD[,3] , col = "antiquewhite4", lty = 1, lwd = 2)
lines(ONSD[,1],ONSD[,4], col = "blue", lty = 2, lwd = 2)
lines(ONSD[,1],ONSD[,5], col = "brown", lty = 3, lwd = 2)
lines(ONSD[,1],ONSD[,6], col = "chartreuse4", lty = 4, lwd = 2)
lines(ONSD[,1],ONSD[,7], col = "darkorchid", lty = 5, lwd = 2)

legend("bottomright", c("tau = 0.05","tau = 0.25", "tau = 0.50", "tau = 0.75","tau = 0.95"),lty =1:5, col = c("antiquewhite4", "blue", "brown","chartreuse4","darkorchid"))

P25=predict(model.q0.25,interval = 'confidence')
View(P25)

names(model.q0.25)



############################

