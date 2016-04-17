baylinreg <- function(Y,X,m_0,S_0,a_0,b_0) {
#
#	Bayesian linear regression model (conjugate prior)
#
# 	INPUTS:
# 	Y					-	Tx1 matrix of system outputs
# 	X					-	Txp-vector of system inputs
# 	m_0				-	px1 matrix, with prior mean 
# 	S_0				-	pxp matrix, with prior precision 
# 	a_0				-	prior shape 
# 	b_0				-	prior scale 
#
#	Author: 		Juan~Carlos Martínez-Ovando
#	Email:		juan.martinez.ovando@itam.mx
#					jc.martinez.ovando@gmail.com
#

#	Initialization
T <- dim(Y)[1]
p <- dim(X)[2]

#	Updating
S_1 <- S_0 + t(X)%*%X
m_1 <- solve(S_1)%*%(S_0%*%m_0 + t(X)%*%Y)
a_1 <- a_0 + T/2
b_1 <- b_0 + 1/2*t(Y-X%*%m_1)%*%Y + 1/2*t(m_0-m_1)%*%S_0%*%m_0

#	Output
output <- list(m_1,S_1,a_1,b_1)
return(output)
}
#
#	END of "baylinreg.R"

bayhierlinreg <- function(hdata,m_0,S_0,a_0,b_0,M) {
#
#	Bayesian linear regression model (conjugate prior)
#
# 	INPUTS:
# 	hdata			-	Tx(1+1+p) data frame of system inputs (a column with consecutive categories for the group id)
# 	m_0				-	px1 matrix, with prior mean 
# 	S_0				-	pxp matrix, with prior precision 
# 	a_0				-	prior shape 
# 	b_0				-	prior scale 
# 	M					-	number of MCMC simulations
#
#	Author: 		Juan~Carlos Martínez-Ovando
#	Email:			juan.martinez.ovando@itam.mx
#						jc.martinez.ovando@gmail.com
#

#	Y.- First column in "hdata"
#	group.- Second colulmn in "hdata"
#	X.- Rest of the columns in "hdata"

#	Initialization
p <- dim(hdata)[2]-2
X <- hdata[,(2+1):(2+p)]
G <- hdata[,2]
Y <- hdata[,1]

J <- max(unique(G))

#	Repositorio
mjh_sim <- array(NaN,c(M,p,J))
mh_sim <- array(NaN,c(M,p))
lambdah_sim <- matrix(NaN,M,1)

#	------------------------------------------------------------
#		Gibbs sampler
#	------------------------------------------------------------
lambda_aux <- rgamma(1, shape=a_0, scale = 1/b_0)
m_aux <- t(as.matrix(mvrnorm(n = 1, mu=m_0, Sigma=solve(lambda_aux*S_0))))
mj_aux <- as.matrix(mvrnorm(n = J, mu=m_0, Sigma=solve(lambda_aux*S_0)))

m <- 1
j <- 1
for(m in 1:M){
	sum_aux <- 0
	#	Simulate latent regression coefficients
	for(j in 1:J){
		Xj <- as.matrix(X[which(hdata[,"clave"]==j),])
		Yj <- as.matrix(Y[which(hdata[,"clave"]==j)])
		#	Updating
		Sj_1 <- S_0 + t(Xj)%*%Xj
		mj_1 <- solve(Sj_1)%*%(S_0%*%m_0 + t(Xj)%*%Yj)
		sum_aux <- sum_aux + 1/2*t(Yj-Xj%*%mj_1)%*%Yj 
		#	Simulation
		mj_aux[j,] <- t(as.matrix(mvrnorm(n = 1, mu=mj_1, Sigma=solve(lambda_aux*Sj_1))))
	  }
	#	Simulate regression coefficients parameters
	#	Updating
	S_1 <- S_0 + (t(mj_aux - t(matrix(colMeans(mj_aux),p,J)))%*%(mj_aux - t(matrix(colMeans(mj_aux),p,J))) )
	m_1 <-  solve(S_1)%*%(S_0%*%m_0 + (t(mj_aux - t(matrix(colMeans(mj_aux),p,J)))%*%(mj_aux - t(matrix(colMeans(mj_aux),p,J))) )%*%as.matrix(colMeans(mj_aux)))
	m_aux <- t(as.matrix(mvrnorm(n = 1, mu=m_1, Sigma=solve(lambda_aux*S_1))))
	#	Updating
	a_1 <- a_0 + T/2
	b_1 <- b_0 + sum_aux + 1/2*t(m_0-m_1)%*%S_0%*%m_0
	#	Simulate scale parameters
	lambda_aux <- rgamma(1, shape=a_1, scale = 1/b_1)
	#  Store samples
	mh_sim[m, ] <- t(m_aux)
	mjh_sim[m, , ] <- t(mj_aux)
	lambdah_sim[m, ] <- lambda_aux 
  }

#  Point estimates	
mjh_mean <- matrix(NaN,p,J)
mh_mean <- matrix(NaN,p,1)
j <- 1
for(j in 1:J){
	mjh_mean[ ,j] <- as.matrix(colMeans(mjh_sim[,,j]))
  }
mh_mean <- as.matrix(colMeans(mh_sim))
lambdah_mean <- mean(lambdah_sim)

#	Output
output <- list(mh_mean,mjh_mean,lambdah_mean,mjh_sim,mh_sim,lambdah_sim)
return(output)
}
#
#	END of "bayhierlinreg.R"
