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
#	Author: 		Juan~Carlos Mart?nez-Ovando
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