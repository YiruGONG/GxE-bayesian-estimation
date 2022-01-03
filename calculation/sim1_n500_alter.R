library(snowfall)
library(coda)
library(MASS)

source('D:/Jiang Lab/GxE/calculation/SIM_regression_alter.R')

Sys.time()
# data generation/simulation
n = 100
p = 4
# beta = c(0.7, 0.1, 0.7, 0.1)
beta = c(0.5,0.5,0.5,0.5)
G = rnorm(n, 0, 1)
E = rnorm(n, 0, 1)
x = matrix(0, ncol=p, nrow=n)
x[,1] = G
x[,2] = G*G-1
x[,3] = E
x[,4] = G*E

g = function(beta, x)
{
  t = x%*%beta
  g = t^2
  return(g)
}
var = 1
y = g(beta, x) + rnorm(n, 0, sqrt(var)) #n*1
y = as.vector(y)

y1 = scale(y)
y1 = as.vector(y1)

sd_y = sd(y)

#-------------------------------------------------------------------------------------
# missing E
E0 = G+rnorm(n, 0, 0.5)        ###???
x1 = matrix(0, ncol=p, nrow=n)
x1[,1] = G
x1[,2] = G*G-1
x1[,3] = E0
x1[,4] = G*E0

data1 = data.frame(y=y1,x=x1)

# initialization
# beta
lm.fit = lm(y~. , data=data1) # [,1:3])
beta0 = lm.fit$coefficients[-1]
beta0 = beta0/(sqrt(sum(beta0^2)))
# beta0 = c(beta0,0,0)
# variance
var0 = sum((lm.fit$residuals)^2)/(n-p-1)
# g
g0 = y1 - rnorm(n, 0, sqrt(var0))           ## y1-rnorm ???

# ESTIMATION
#----------------------------------------------------------------------------------------------------------
# Method2: proposal=='t_proj'
# parameter
var_b = 1*10^(-5)
var_E = 7.5*10^(-4)
var_tau = 2
var_l = 0.2
sigma_b0 = rep(1, p)
a = 0.5
b = 0.5
tau0 = 5
l0 = 0.1
ite = 2000

result = mh_gib1(y1,G,E0,g0,beta0,var0,tau0,l0,var_b=var_b,sigma_b0=sigma_b0,a=a,b=b,var_E=var_E,var_tau=var_tau,var_l=var_l,ite=ite,proposal='vmf')
result[[7]]
beta = result[[2]]

Sys.time()

pdf(file='D:/Jiang Lab/GxE/results/calculation_vmf_minb_orig.pdf',width=10,height=10)
par(mfrow=c(4,1))
plot(beta[,1], type = "l",col = "blue", xlab = "Iteration", ylab = "Estimation", main=expression(beta[1]))
plot(beta[,2], type = "l",col = "blue", xlab = "Iteration", ylab = "Estimation", main=expression(beta[2]))
plot(beta[,3], type = "l",col = "blue", xlab = "Iteration", ylab = "Estimation", main=expression(beta[3]))
plot(beta[,4], type = "l",col = "blue", xlab = "Iteration", ylab = "Estimation", main=expression(beta[4]))
dev.off()

##other observations
Emean = apply(result[[1]],1,mean)
Esd = apply(result[[1]],1,var)
par(mfrow=c(2,1))
plot(Emean,type='l',col='blue',xlab='Iteration',ylab='Estimation',main='E mean')
plot(Esd,type='l',col='blue',xlab='Iteration',ylab='Estimation',main='E sd')

