#
#	EST-46114:	Inferencia Bayesiana en Alta Dimension (Maestria en Ciencia de Datos)
#	Autor: 		Juan~Carlos Martinez Ovando
#	Email:		juan.martinez.ovando@itam.mx
#				jc.martinez.ovando@gmail.com
#	
#	Seleccion Estocastica de Variables
#	== Modelos Lineales Generalizados: Probit ==
#

rm(list = ls())
rutawork = 'C:/JCMO.Academia/@Cursos/2016-I_Inferencia Bayesiana en Alta Dimension/_data&code/SSVS/'
#rutawork <- '/media/jcmo/ROCADATA/JCMO.Academia/@Cursos/2016-I_Inferencia\ Bayesiana\ en\ Alta Dimension/_data&code/SSVS/'

#	----------------------------------------------------
#		Libraries
#	----------------------------------------------------
install.packages('mvtnorm')
library(mvtnorm)

#	----------------------------------------------------
#		Datos
#	----------------------------------------------------
dde <- read.csv(paste(rutawork,'EST46114_SSVS_Data.csv',sep = ""), header = TRUE, sep = ",", quote="\"", dec=".", fill = TRUE)

# 	Keys
# x = dde 	# dosis
# Y = 0/1 	# variable indicadora sobre embarazo premanuturoindicato
# z1-z5 	# covariables (variables confundidas)

# 	Normalizacion 
dde$x<- (dde$x - mean(dde$x))/sqrt(var(dde$x)) 

#	--------------------------------------------			
#		Analisis frecuentista
#	--------------------------------------------			
X <- cbind(1,dde$x,dde$z1,dde$z2,dde$z3,dde$z4,dde$z5)
Y <- dde$y
dde_mle<- glm(Y ~ -1+X, family=binomial("probit"))

# Summary table 
#Coefficients:
#            		Estimate	Std. Error 		z value		Pr(>|z|)
#(Intercept)	-1.08068	0.04355 		-24.816  	< 2e-16 ***
#dde			 0.17536	0.02909   		  6.028 		1.67e-09 ***
#z1				-0.12817	0.03528		 -3.633 	0.000280 ***
#z2				 0.11097	0.03366		  3.297 		0.000978 ***
#z3				-0.01705	0.03405		 -0.501 	0.616659
#z4				-0.08216	0.03576		 -2.298 	0.021571 *
#z5				 0.05462 	0.06473		  0.844 		0.398721

#	Estimador maximo verosimil
beta_mle<- dde_mle$coef

#	--------------------------------------------			
#		Analisis bayesiano
#	--------------------------------------------			

# 		Especificacion inicial 
beta0 <- rep(0,7)
Pbeta0 <- 0.25*diag(7)

#		MCMC
beta <- rep(0,7)	# Valor inicial de la cadena de Markov
n <- nrow(dde)		# Numero de observaciones/individuos
z <- rep(0,n)		# Valores iniciales de las variables latentes (modelo probit)
G <- 1000			# Numero de iteraciones del MCMC

# 		Gibbs sampler
gt <- 1
for(gt in 1:G){
  eta<- X%*%beta 	# Predictor lineal (variable auxiliar)
  
  #	Muestreo de la distribucion truncada para Z 
  #	de las distribuciones condicionales completas
  z[Y==0] <- qnorm(runif(sum(1-Y), 0,pnorm(0,eta[Y==0],1)), eta[Y==0],1)
  z[Y==1] <- qnorm(runif(sum(Y), pnorm(0,eta[Y==1],1),1), eta[Y==1],1)

  #	Muestreo de la distribucion final completa para beta
  Vbeta <- solve(Pbeta0 + t(X)%*%X)
  Ebeta <- Vbeta%*%(Pbeta0%*%beta0 + t(X)%*%z)
  beta <- c(rmvnorm(1,Ebeta,Vbeta))

  #	Output
  write(t(beta),file=paste(rutawork,'EST46114_SSVS_beta.out',sep = ""), ncol=7, append=T)
  print(c(gt,round(beta*100)/100))
}

#	Trace-plots para beta[1] y beta[2]
pdf(paste(rutawork,'EST46114_SSVS_beta_traces.pdf',sep = ""),width=7,height=5)
par(mfrow=c(2,1))
beta_out<- matrix(scan(paste(rutawork,'EST46114_SSVS_beta.out',sep = "")), ncol=7, byrow=T)
plot(beta_out[,1],type="l",xlab="iteration",
		ylab="intercept (beta_1)")
plot(beta_out[,2],type="l",xlab="iteration",ylab="slope (beta_2)")
dev.off()

#	Distribucion final marginal para el coeficiente de pendiente de DDE
pdf(paste(rutawork,'EST46114_SSVS_beta_marginal_slope.pdf',sep = ""),width=7,height=5)
slp <- beta_out[101:1000,2]
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
plot(density(slp),type="l",xlab="DDE slope (beta_1)",ylab="Posterior Density",
     cex=1.2)
abline(v=mean(slp))
abline(v=0.17536139,lwd=2.5) # MLE
abline(v=0.17536139 + 1.96*0.02909*c(-1,1),lwd=2.5,lty=2)
abline(v=quantile(slp,probs=c(0.025,0.975)),lty=2)
dev.off()

#	Curva de respuesta a dosis DDE para nacimiento prematuro
pdf(paste(rutawork,'EST46114_SSVS_beta_doseresponse.pdf',sep = ""),width=7,height=5)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
xg<- seq(-2,2,length=100)  	# Grid
beta<- beta_out[101:1000,] 	# Datos sin periodo de calentamiento
post<- matrix(0,100,4)
for(i in 1:100){
  post[i,1]<- mean(pnorm(beta[,1] + xg[i]*beta[,2]))
  post[i,2:3]<- quantile(pnorm(beta[,1] + xg[i]*beta[,2]),probs=c(0.025,0.975))
  post[i,4]<- pnorm(-1.08068 + xg[i]*0.17536139)
}
xtrue<- xg*sqrt(var(dde$x)) + mean(dde$x)
plot(xtrue,post[,1],xlab="Serum DDE (mg/L)",ylab="Pr(Preterm birth)",cex=1.2,
     ylim=c(0,max(post)), type="l")
lines(xtrue,post[,2],lty=2)
lines(xtrue,post[,3],lty=2)
lines(xtrue,post[,4],lwd=2.5)
dev.off()

#	Resumen de la distrinbucion final de los coeficientes de regresion
table1<- matrix(0,7,5)
for(i in 1:7){
  table1[i,]<- c(mean(beta[,i]),median(beta[,i]),sqrt(var(beta[,i])),
                 quantile(beta[,i],probs=c(0.025,0.975)))
}
table1<- round(table1*100)/100
write.csv(t(table1),file=paste(rutawork,'EST46114_SSVS_beta_summary.csv',sep = ""))

#	---------------------------------------------------------------------			
#		Selección estocástica de variables via el Gibbs sampler
#	---------------------------------------------------------------------			

ddemle <- glm(Y ~ -1+X, family=binomial("probit"))
betamle <- ddemle$coef 	#	Estimación por máxima verosimilitud (solo para el valor incial del MCMC)

#	Especificacion inicial
p <- ncol(X)        	# Numero de covariables (predictores)
p0 <- rep(0.5,p)    	# Distribucion inicial para los predictores excluidos
b0 <- rep(0,p)      	# Hiperparametro para el predictor seleccionado
s0 <- rep(2,p)      	# Varianza del componente normal latente
	
#	Stochastic search variable selection via data-augmentation
beta <- betamle
n <- nrow(dde)   	# Numero de individuos
z <- rep(0,n)    	# Valores iniciales de variables latentes
G <- 5000        	# Numero de iteraciones MCMC

# -------------------------------------
#	Stochastic search Gibbs sampler
# -------------------------------------
gt <- 1
for(gt in 1:G){
  #		Data augmentation step
  eta<- X%*%beta 	# Predictor lineal
  
  # Muestreo de la distribucion normal truncada para latentes
  # de la distribucion condicional completa
  z[Y==0]<- qnorm(runif(sum(1-Y),0,pnorm(0,eta[Y==0],1)),eta[Y==0],1)
  z[Y==1]<- qnorm(runif(sum(Y),pnorm(0,eta[Y==1],1),1),eta[Y==1],1)

  # Actualizacion de coeficientes de regresion (uno a la vez, en lugar de por bloque como en el algoritmo anterior)
  j <-1
  for(j in 1:p){
    #	Varianza final condicional para beta_j bajo la distribucion inicial normal
    Vj<- 1/(s0[j]^{-2} + sum(X[,j]^2))
    Ej<- Vj*sum(X[,j]*(z-X[,-j]%*%beta[-j]))
    
	#	Probabilidad condicional de incluir el predictor 'j'
    phat<- 1/(1 + p0[j]/(1-p0[j]) * dnorm(0,Ej,sqrt(Vj))/dnorm(0,b0[j],s0[j]) )          
    m<- rbinom(1,1,phat)
    beta[j]<- m*rnorm(1,Ej,sqrt(Vj))
  }
  
  #	Output
  write(t(beta),file=paste(rutawork,'EST46114_SSVS_beta_ss.out',sep = ""), ncol=7, append=T)
  print(c(gt,round(beta*100)/100))
}

#	Trace-plots para beta[1] y beta[2]
beta<- matrix(scan(paste(rutawork,'EST46114_SSVS_beta_ss.out',sep = "")), ncol=7, byrow=T)
# 2^p = Numero de posibles modelos
# En este caso, 2^p = 2^7 = 128

pdf(paste(rutawork,'EST46114_SSVS_beta_traces_ss.pdf',sep = ""),width=7,height=5)

# 	Grafica de las iteraciones de Gibbs
par(mfrow=c(3,2))
ylb=c("intercept","dde","z1","z2","z3","z4","z5")
for(j in 2:7){
  print(j)
  plot(beta[,j],xlab="iteration",ylab=ylb[j])
}
dev.off()

#	Grafica de la curva de respuesta a DDE para nacimiento prematuro
pdf(paste(rutawork,'EST46114_SSVS_beta_doseresponse_ss.pdf',sep = ""),width=7,height=5)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
xg<- seq(-2,2,length=100)  	# Grid
beta<- beta[1001:nrow(beta),] 	# Periodo de calentamiento
post<- matrix(0,100,4)
for(i in 1:100){
  post[i,1]<- mean(pnorm(beta[,1] + xg[i]*beta[,2]))
  post[i,2:3]<- quantile(pnorm(beta[,1] + xg[i]*beta[,2]),probs=c(0.025,0.975))
  post[i,4]<- pnorm(-1.08068 + xg[i]*0.17536139)
}
xtrue<- xg*sqrt(var(dde$x)) + mean(dde$x)
plot(xtrue,post[,1],xlab="Serum DDE (mg/L)",ylab="Pr(Preterm birth)",cex=1.2,
     ylim=c(0,max(post)), type="l")
lines(xtrue,post[,2],lty=2)
lines(xtrue,post[,3],lty=2)
lines(xtrue,post[,4],lwd=2.5)
dev.off()

#	Resumen de coeficientes de regresion
#	Incluye:
#		- media posterior
#		- mediana
#		- desviacion estandar
#		- 95% intervalo de confianza
#		- Pr(beta_j=0|datos)
table1<- matrix(0,7,6)
for(i in 1:7){
  table1[i,]<- c(mean(beta[,i]),median(beta[,i]),sqrt(var(beta[,i])),
                 quantile(beta[,i],probs=c(0.025,0.975)),
		 length(beta[beta[,i]==0,i])/nrow(beta))
}
table1<- round(table1*100)/100
write(t(table1),file=paste(rutawork,'EST46114_SSVS_beta_summary_ss.csv',sep = ""),ncol=6)

#	Probabilidades finales para los mejores modelos visitados por el algoritmo
Mout <- beta 
Mout <- matrix(as.numeric(I(Mout==0)),nrow(Mout),ncol(Mout)) 
Mindex <- Mout[1,] 		   # Different models sampled starting with first returns 1 if m1 and m2 are the same model
ind<- function(m1,m2){
  as.numeric(all(I(m1==m2)))
}
Im<- apply(Mout,1,ind,m2=Mout[1,]) # Indicadora para las muestras del primer modelo
Nm<- sum(Im)                       # Numero de muestral del primer modelo
Mout<- Mout[Im==0,]
repeat{
  if(length(Mout)==7) Mout<- matrix(Mout,1,7)
  Mindex<- rbind(Mindex,Mout[1,])
  Im<- apply(Mout,1,ind,m2=Mout[1,])
  Nm<- c(Nm,sum(Im))
  print(sum(Nm))
  if(sum(Nm)==nrow(beta)){ 
    break
  } else Mout<- Mout[Im==0,]
}
#	Ordenando los modelos visitados en terminos de sus probabilidades (ordendecreciente)
Pm <- Nm/sum(Nm)
ord <- order(Pm)
Pm <- Pm[rev(ord)]
Mindex <- Mindex[rev(ord),]
table2 <- cbind(Pm,Mindex)
write(t(table2),file=paste(rutawork,'EST46114_SSVS_beta_pmodel_ss.csv',sep = ""),ncol=8)

#
#	-- FIN: EST46114_SSVS_Script.r --