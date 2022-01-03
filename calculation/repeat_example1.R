setwd("D:/Jiang Lab/GxE")

source('./calculation/repeat_regression.R')

### example 1
################ data generation ###################
# n = 100
# x = cbind(runif(n,-3,5),rnorm(n,0,9))
# beta = c(0.5,0.5*sqrt(3))
# sig = 0.01
# gf = function(t) {2*t^3 + 3*sin((0.5*t)^3)}
# g = gf(x%*%beta)
# y = g + rnorm(n,0,sqrt(sig))

##try sample 2
n = 100
p = 4
beta = c(0.8,0.4,-0.4,-0.2)
x = matrix(0,ncol=p,nrow=n)
x[,1]=runif(n,-3,5)
x[,2]=rnorm(n,0,9)
x[,3]=rchisq(n,3)
x[,4]=x[,2]*x[,3]
x = apply(x,2,scale)

g = function(beta,x){
  t= x%*%beta
  g = 5*cos(t)+exp(-t^2)
  return(g)
}
sig = 0.1
y = g(beta,x) + rnorm(n,0,sqrt(var))
y = as.vector(y)

y1 = scale(y)
y1 = as.vector(y1)
sd_y = sd(y)

################## estimation ######################
n = nrow(x)
p = ncol(x)

data = data.frame(y=y1,x=x)
# colnames(data) = c('y','x1','x2')
model = lm(y~.,data=data)
beta0 = model$coefficients[-1]
beta0 = beta0/(sqrt(sum(beta0^2)))  ##standardized
sig0 = sum((model$residuals)^2)/(n-p-1)
g0 = y - rnorm(n,0,sqrt(sig0))
# g0 = x%*%beta0
# sig0 = sum((y-g0)^2)/(n-p-1)

a=0.5
b=0.5
var_tau = 2
var_l = 0.2
tau0 = 5
l0 = 0.1
ite = 1000

results = estimate_bayes(y,x,beta0,tau0,l0,g0,sig0,var_tau,var_l,a,b,iteration=ite,method='vmf')
results[[6]]
beta = results$beta

pdf(file='D:/Jiang Lab/GxE/results/calculation_vmf_minb_orig.pdf',width=10,height=10)
par(mfrow=c(2,1))
plot(beta[,1], type = "l",col = "blue", xlab = "Iteration", ylab = "Estimation", main=paste(expression(beta[1],method1)))
plot(beta[,2], type = "l",col = "blue", xlab = "Iteration", ylab = "Estimation", main=paste(expression(beta[2],method1)))
dev.off()

