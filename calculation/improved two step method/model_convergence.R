setwd('D:/Jiang Lab/GxE')

library( PearsonDS )
# library(RNOmni)
library(doParallel)
library(foreach)

library(reshape2)
library(gridExtra)
library(ggplot2)
# library(gsl)

source('./calculation/two step/estimate_basic.R')


####### estimation for diff Es #########
m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0
b1   = sqrt(.3)  # beta
c1   = sqrt(.05)  # gamma
sig  = 0.5
# pow  = 1
# nEs   = c(seq(2,30,2),seq(35,50,5))
nEs  = c(1,seq(5,50,5))
# nE = 1

cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)
# clusterExport(cl,c(),envir = environment())
Sys.time()

convergence_Ent_MLE = foreach(nE=nEs) %dopar% {
  # start=Sys.time()
  # alld = lapply(t(as.matrix(d)),function(d1){
  reps = sapply(1:500,function(rep){
    maf  =  matrix( runif( m ), nrow = 1 )
    G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    eff  =  matrix( runif( m ), ncol = 1 )
    eff  =  eff / sqrt( sum( eff^2 ) )
    GRS  =  scale( G %*% eff )
    
    # idxE = sample(1:8,nE,replace=T)
    # idxE = rep(1:8,5)[1:nE]
    # Es   = sapply(idxE,function(x) {
    #   if (x==1) {runif(n,-sqrt(3),sqrt(3))
    #   }else if (x==2) { rnorm(n)
    #   }else if (x==3) { (rchisq(n,3)-3)/sqrt(6)
    #   }else if (x==4) { rt(n,3,ncp=0)/sqrt(3)
    #   }else if (x==5) { rbinom(n,4,0.5)-2
    #   }else if (x==6) { rpois(n,1)-1
    #   }else if (x==7) { scale(rbeta(n,0.5,0.5))
    #   }else if (x==8) { scale(rgamma(n,0.5)) }
    #   })
    Es = sapply(1:nE,function(x) scale(rt(n,3,ncp=0)/sqrt(3)))
    # Es = apply(Es,2,scale)
    coeff.E = rep(0.5,nE)
    E = Es %*% coeff.E / sqrt(sum(coeff.E^2))
    
    noi = rnorm(n,0,1)
    z    =  a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi
    y = z - mean(z)
    
    # Estimate interaction effect for GRS
    pars = estimate_basic( y, GRS,method='MLE' )
    return( pars )
  })
  mean = apply(reps,1,mean)
  var  = apply(reps,1,var)
  return( list(mean=mean,var=var,reps=reps) )
}
stopCluster(cl)
Sys.time()
save(convergence_Ent,file='./data/new/convergence_Ent.Rdata')

#### plots #######
data = convergence_Ent_MLE
data = lapply(convergence,function(x) {
  data2 = list()
  for (n in 1:500) {
    est = x[(3*n-2):(3*n)]
    data2[[n]] = est
  }
  return(data2)
})

alpha1 = sapply(data,function(x) x$reps['alpha1',])
alpha2 = sapply(data,function(x) x$reps['alpha2',])
beta = sapply(data,function(x) x$reps['beta',])
gamma = sapply(data,function(x) x$reps['gamma',])
sigma = sapply(data,function(x) x$reps['sigma',])
colnames(alpha1) <- colnames(alpha2) <- colnames(beta) <- colnames(gamma) <- colnames(sigma) <- nEs


# alpha1 = sapply(data,function(x) sapply(x,function(xx) xx$pars['alpha1']))
# alpha2 = sapply(data,function(x) sapply(x,function(xx) xx$pars['alpha2']))
# beta = sapply(data,function(x) sapply(x,function(xx) xx$pars['beta']))
# gamma = sapply(data,function(x) sapply(x,function(xx) xx$pars['gamma']))
# sigma2 = sapply(data,function(x) sapply(x,function(xx) xx$pars['sigma2']))
# colnames(alpha1) <- colnames(alpha2) <- colnames(beta) <- colnames(gamma) <- colnames(sigma2) <- nEs

pdf('./results/new/convergence_Ent_MLE.pdf',width = 12,height=12)
par(mfrow=c(3,2))
boxplot(alpha1,xlab='different number of En',ylab='Estimation for alpha1',main='alpha1',las=1)
abline(h=a1,lty=2)
boxplot(alpha2,xlab='different number of En',ylab='Estimation for alpha2',main='alpha2',las=1)
abline(h=a2,lty=2)
boxplot(beta,xlab='different number of En',ylab='Estimation for beta',main='beta',las=1)
abline(h=b1,lty=2)
boxplot(gamma,xlab='different number of En',ylab='Estimation for gamma',main='gamma',las=1)
abline(h=c1,lty=2)
boxplot(sigma,xlab='different number of En',ylab='Estimation for sigma',main='sigma',las=1)
abline(h=sig,lty=2)
boxplot(gamma,xlab='different number of En',ylab='Estimation for gamma',main='zoom - gamma',las=1,ylim=c(0.15,0.3))
abline(h=c1,lty=2)
dev.off()


### calculate bias and var for gamma
bias = apply(gamma,2,mean)-c1
var  = apply(gamma,2,var)

pdf('./results/new/convergence_bias_var.pdf',width = 8,height=4)
par(mfrow=c(1,2))
plot(bias~nEs,type='l',xlab='number of En',las=1)
plot(var~nEs,type='l',xlab='number of En',ylab='Variance',las=1)
dev.off()

## view convergence of 1

####### estimation for diff distribution #########
m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0
b1   = sqrt(.3)  # beta
c1   = sqrt(.05)  # gamma
sig  = 0.5
# pow  = 1
# nEs   = c(1:4, seq(5,30,5))
nE = 1

cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)
# clusterExport(cl,c(),envir = environment())
Sys.time()

convergence_nE1 = foreach(idxE=1:8) %dopar% {
  # start=Sys.time()
  # alld = lapply(t(as.matrix(d)),function(d1){
  reps = lapply(1:500,function(rep){
    maf  =  matrix( runif( m ), nrow = 1 )
    G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    eff  =  matrix( runif( m ), ncol = 1 )
    eff  =  eff / sqrt( sum( eff^2 ) )
    GRS  =  scale( G %*% eff )
    
    # idxE = sample(1:8,nE,replace=T)
    # idxE = rep(1:8,5)[1:nE]
    E   = sapply(idxE,function(x) {
      if (x==1) {runif(n,-sqrt(3),sqrt(3))
      }else if (x==2) { rnorm(n)
      }else if (x==3) { (rchisq(n,3)-3)/sqrt(6)
      }else if (x==4) { rt(n,3,ncp=0)/sqrt(3)
      }else if (x==5) { rbinom(n,4,0.5)-2
      }else if (x==6) { rpois(n,1)-1
      }else if (x==7) { scale(rbeta(n,0.5,0.5))
      }else if (x==8) { scale(rgamma(n,0.5)) }
    })
    # Es = apply(Es,2,scale)
    
    noi = rnorm(n,0,1)
    z    =  a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi
    y = z - mean(z)
    
    # Estimate interaction effect for GRS
    pars = estimate_basic( y, GRS,method='OLS' )
    return( list(E=E,pars=pars) )
  })
  return(reps)
}
stopCluster(cl)
Sys.time()
save(convergence_nE1,file='./data/new/convergence_nE1.Rdata')

data = convergence_nE1
alpha1 = sapply(data,function(x) sapply(x,function(xx) xx[['pars']]['alpha1']))
alpha2 = sapply(data,function(x) sapply(x,function(xx) xx[['pars']]['alpha2']))
beta = sapply(data,function(x) sapply(x,function(xx) xx[['pars']]['beta']))
gamma = sapply(data,function(x) sapply(x,function(xx) xx[['pars']]['gamma']))
sigma = sapply(data,function(x) sapply(x,function(xx) xx[['pars']]['sigma']))
colnames(alpha1) <- colnames(alpha2) <- colnames(beta) <- colnames(gamma) <- colnames(sigma) <- c(1:8)

pdf('./results/new/convergence_nE1.pdf',width = 10,height=16)
par(mfrow=c(4,2))
boxplot(alpha1,xlab='different distribution of E',ylab='Estimation for alpha1',main='alpha1',las=1)
abline(h=a1,lty=2)
boxplot(alpha2,xlab='different distribution of E',ylab='Estimation for alpha2',main='alpha2',las=1)
abline(h=a2,lty=2)
boxplot(beta,xlab='different distribution of E',ylab='Estimation for beta',main='beta',las=1)
abline(h=b1,lty=2)
boxplot(gamma,xlab='different distribution of E',ylab='Estimation for gamma',main='gamma',las=1)
abline(h=c1,lty=2)
boxplot(beta,xlab='different distribution of E',ylab='Estimation for beta',main='zoom - beta',las=1,ylim=c(0.4,0.7))
abline(h=b1,lty=2)
boxplot(gamma,xlab='different distribution of E',ylab='Estimation for gamma',main='zoom - gamma',las=1,ylim=c(0.1,0.3))
abline(h=c1,lty=2)
boxplot(sigma,xlab='different distribution of E',ylab='Estimation for sigma',main='sigma',las=1)
abline(h=sig,lty=2)
boxplot(sigma,xlab='different distribution of E',ylab='Estimation for sigma',main='zoom-sigma',las=1,ylim=c(0.4,0.6))
abline(h=sig,lty=2)
dev.off()

### bias,var v.s. e3,e4
bias = apply(gamma,2,mean)-c1
var  = apply(gamma,2,var)
e3s   = sapply(data, function(x) sapply(x,function(xx) mean(xx$E^3)))
e4s   = sapply(data, function(x) sapply(x,function(xx) mean(xx$E^4)))
e3 = apply(e3s,2,mean)
e4 = apply(e4s,2,mean)
all = cbind(bias,var,e3,e4)
cor = cor(all)


par(mfrow=c(2,2))
plot(bias~e3)
plot(bias~e4)
plot(var~e3)
plot(var~e4)

# try1 = sapply(data[[1]],function(x) c(x[[1]],x[[3]]))
# try1 = t(try1)
# colnames(try1)[1] = 'distribution'
# # try1 = as.data.frame(try1)
# # try1$distribution = as.factor(try1$distribution)
# 
# parname = c('alpha1','alpha2','beta','gamma','sigma')
# parvalue = c(a1,a2,b1,c1,sig)
# names(parvalue) = parname
# 
# plots=list()
# for (idxp in 1:5){
#   par = parname[idxp]
#   try2 = cbind(distribution=try1[,1],value=try1[,(idxp+1)])
#   try2 = as.data.frame(try2)
#   try2$distribution = as.factor(try2$distribution)
#   plots[[idxp]] <- ggplot(try2,aes(x=distribution,y=value))+geom_boxplot()+
#     geom_hline(yintercept = parvalue[par],linetype='dashed',color='red')+
#     labs(x='different distribution of E',y=par ,title='gamma^2=0.025')
# }
# 
# grid.arrange(grobs=plots,nrow=3)
# 
# 
# # bias2 = apply(alpha1,2,mean)-a1
# # var2 = apply(alpha1,2,'var')
# # rbind(bias1,bias2)
# bias11 = aggregate(gamma~distribution,data=try1,mean)$gamma-c1
# var11 = aggregate(gamma~distribution,data=try1,'var')$gamma
# e3 = apply(Es^3,2,mean)
# e4 = apply(Es^4,2,mean)
# test = cbind(bias11,var11,e3,e4)
# plot(bias11~e3)
# plot(bias11~e4)
# plot(var11~e3)
# plot(var11~e4)
