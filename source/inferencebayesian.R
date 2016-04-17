library(sqldf)
library(rjags)
library(R2jags)
library(ggplot2)
# ESTIMATION OF A MEAN
data <- read.csv("~/Desktop/Fun_Estadistica/EST4611_AbortoCrimen.csv")
#milk<-read.table("http://allman.rhon.itam.mx/~lnieto/index_archivos/milk.txt",header=TRUE)
data$t<-66:99
n<-nrow(data)
or<-order(data$xxprison)
plot(data$xxprison[or],data$lpc_murd[or],type="l")
text(data$xxprison[or],data$lpc_murd[or],labels=data$t[or],cex=0.5,col=2)
plot(data$t,data$lpc_murd,type="l")
plot(data$t,data$xxprison,type="l")


#-Defining data-
m<-2
data<-list("n"=n,"y"=data$lpc_murd,"x"=data$xxprision,"t"=data$t) #"m"=m,
#data<-list("n"=n,"y"=c(milk$y[1:(n-2)],NA,NA),"x"=milk$x/max(milk$x),"t"=milk$t/max(milk$t))

#-Defining inits-
inits<-function(){list(beta=0,tau=1,yf1=rep(1,n))}
#inits<-function(){list(beta=rep(0,5),tau=1,yf1=rep(1,n))}
#inits<-function(){list(beta=rep(0,n+m),tau.y=1,tau.b=1,yf1=rep(0,n+m))}
#inits<-function(){list(beta=rep(0,n+m),tau.y=1,yf1=rep(0,n+m))}
#inits<-function(){list(beta=rep(0,n+m),tau.y=1,tau.b=1,yf1=rep(0,n+m),g=0)}
#inits<-function(){list(beta=rep(0,n),tau.y=1,tau.b=1,yf1=rep(0,n),g=1)}

#-Selecting parameters to monitor-
parameters<-c("beta","tau","yf1")
#parameters<-c("beta","tau.y","tau.b","yf1","g")
#parameters<-c("beta","tau.y","tau.b","yf1")

#-Running code-
#JAGS
model <- function(){
  #Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(mu[i],tau)
    mu[i]<-beta[1]+beta[2]*x[i]+beta[4]*t[i]
  }
  #Priors 
  for (j in 1:3) { beta[j] ~ dnorm(0,0.001) }
  tau ~ dgamma(0.001,0.001)
  
  #Prediction
for (i in 1:n) { yf1[i] ~ dnorm(mu[i],tau) }
}


sim<-jags(data,inits,parameters,model,#.file="Ej7a.txt",
               n.iter=50000,n.chains=1,n.burnin=5000)

#-Monitoring chain-

#Traza de la cadena
traceplot(ej6.sim)

#Cadena

#OpenBUGS
out<-ej7c.sim$sims.list

#JAGS
out<-ej7a.sim$BUGSoutput$sims.list

z<-out$tau
par(mfrow=c(2,2))
plot(z,type="l")
plot(cumsum(z)/(1:length(z)),type="l")
hist(z,freq=FALSE)
acf(z)

#Resumen (estimadores)
#OpenBUGS
out.sum<-ej7c.sim$summary

#JAGS
out.sum<-ej7a.sim$BUGSoutput$summary

print(out.sum)
head(out.sum)

#Probabilidades
z<-out$beta[,1]
prob(z)

#DIC
out.dic<-ej7c.sim$DIC
out.dic<-ej7a.sim$BUGSoutput$DIC
print(out.dic)

#Predictions
out.yf<-out.sum[grep("yf1",rownames(out.sum)),]
ymin<-min(milk$y,out.yf[,c(1,3,7)])
ymax<-max(milk$y,out.yf[,c(1,3,7)])
xmin<-min(milk$t)
xmax<-max(milk$t+m)

#x vs. y
par(mfrow=c(1,1))
plot(milk$x,milk$y,type="p",col="grey50",ylim=c(ymin,ymax))
points(milk$x,out.yf[,1],col=2,pch=16,cex=0.5)
segments(milk$x,out.yf[,3],milk$x,out.yf[,7],col=2)

#t vs y
par(mfrow=c(1,1))
plot(milk$t,milk$y,type="b",col="grey80",ylim=c(ymin,ymax),xlim=c(xmin,xmax))
lines(milk$t,out.yf[1:n,1],col=2)
lines(milk$t,out.yf[1:n,3],col=2,lty=2)
lines(milk$t,out.yf[1:n,7],col=2,lty=2)
lines(milk$t[n]:(milk$t[n]+m),out.yf[n:(n+m),1],col=4)
lines(milk$t[n]:(milk$t[n]+m),out.yf[n:(n+m),3],col=4,lty=2)
lines(milk$t[n]:(milk$t[n]+m),out.yf[n:(n+m),7],col=4,lty=2)

#betas
out.beta<-out.sum[grep("beta",rownames(out.sum)),]
ymin<-min(out.beta[,c(1,3,7)])
ymax<-max(out.beta[,c(1,3,7)])
plot(out.beta[,1],type="l",ylim=c(ymin,ymax))
lines(out.beta[,3],lty=2)
lines(out.beta[,7],lty=2)