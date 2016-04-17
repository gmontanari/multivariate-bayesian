#data(package = .packages(all.available = TRUE))

library(hglm)
vignette('hglm')
data(semiconductor)

m11 <- hglm(fixed = y ~ x1 + x3 + x5 + x6,
            random = ~ 1|Device,
            family = Gamma(link = log),
            disp = ~ x2 + x3, data = semiconductor)
summary(m11)
plot(m11, cex = .6, pch = 1,
     cex.axis = 1/.6, cex.lab = 1/.6,
     cex.main = 1/.6, mar = c(3, 4.5, 0, 1.5))
# ------------------- #
# redo it using hglm2 #
# ------------------- #
m12 <- hglm2(y ~ x1 + x3 + x5 + x6 + (1|Device),
             family = Gamma(link = log),
             disp = ~ x2 + x3, data = semiconductor)
summary(m12)
# -------------------------- #
# redo it using matrix input #
# -------------------------- #
attach(semiconductor)
m13 <- hglm(y = y, X = model.matrix(~ x1 + x3 + x5 + x6),
            Z = kronecker(diag(16), rep(1, 4)),
            X.disp = model.matrix(~ x2 + x3),
            family = Gamma(link = log))
summary(m13)
# --------------------- #
# verbose & likelihoods #
# --------------------- #
m14 <- hglm(fixed = y ~ x1 + x3 + x5 + x6,
            random = ~ 1|Device,
            family = Gamma(link = log),
            disp = ~ x2 + x3, data = semiconductor,
            verbose = TRUE, calc.like = TRUE)
summary(m14)
# --------------------------------------------- #
# simulated example with 2 random effects terms #
# --------------------------------------------- #
## Not run:
set.seed(911)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
z1 <- factor(rep(LETTERS[1:10], rep(10, 10)))
z2 <- factor(rep(letters[1:5], rep(20, 5)))
Z1 <- model.matrix(~ 0 + z1)
Z2 <- model.matrix(~ 0 + z2)
u1 <- rnorm(10, 0, sqrt(2))
u2 <- rnorm(5, 0, sqrt(3))
y <- 1 + 2*x1 + 3*x2 + Z1%*%u1 + Z2%*%u2 + rnorm(100, 0, sqrt(exp(x3)))
dd <- data.frame(x1 = x1, x2 = x2, x3 = x3, z1 = z1, z2 = z2, y = y)
(m21 <- hglm(X = cbind(rep(1, 100), x1, x2), y = y, Z = cbind(Z1, Z2),
             RandC = c(10, 5)))
summary(m21)
plot(m21)

m31 <- hglm2(y ~ x1 + x2 + (1|z1) + (1|z2), disp = ~ x3, data = dd)
print (m31)
summary(m31)
plot(m31)

#likelihiood ratio test
m20 <- hglm(X = cbind(rep(1, 100), x1, x2), y = y, Z = Z1,
            calc.like = TRUE)
lrt(m20)
m21 <- hglm(X = cbind(rep(1, 100), x1, x2), y = y, Z = cbind(Z1, Z2),
            RandC = c(10, 5), calc.like = TRUE)
lrt(m20, m21)


##EJEMPLOS
library(MASS)
data(bacteria)
g1 <- hglm(fixed = y ~ week, random = ~1 | ID, data = bacteria, family = binomial(link = logit))
summary(g1)
plot(g1, device = "pdf", name = "fig_bacteria")


##linear mixed model
n.cluster <- 5
n.per.cluster <- 20
sigma2_u <- 0.2
sigma2_e <- 1
beta.disp <- 1
mu <- 0
n <- n.cluster * n.per.cluster
set.seed(1234)
X <- matrix(1, n, 1)
Z <- diag(n.cluster) %x% rep(1, n.per.cluster)
a <- rnorm(5, 0, sqrt(sigma2_u))
X_d <- matrix(1, n, 2)
X_d[, 2] <- rbinom(n, 1, 0.5)
e <- rnorm(n, 0, sqrt(sigma2_e * exp(beta.disp * X_d[, 2])))
y <- mu + Z %*% a + e
simul1 <- hglm(y = y, X = X, Z = Z, X.disp = X_d)
summary(simul1)





