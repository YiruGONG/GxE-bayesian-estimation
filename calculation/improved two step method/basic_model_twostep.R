setwd('D:/Jiang Lab/GxE')

library( PearsonDS )
# library(RNOmni)
library(doParallel)
library(foreach)

library(reshape2)
library(gridExtra)
library(ggplot2)
# library(gsl)

# source('./calculation/two step/step1_fit.R')
# source('./calculation/two step/step2_fit.R')
# source('./calculation/two step/OLS.R')
source('./calculation/two step/estimate_basic.R')

.single_gxe  =  function( y,
                          grs,
                          params ){
  par1  =  optim( params[1:2],
                    step1_fit,
                    y = y, grs = grs )$par
  par2  = optim(  params[3:5],
         step2_fit,
         y = y, grs = grs, a = par1 )$par
  c(par1, par2)
}

estimate_basic = function(y,grs,method){
  if (method=='MLE'){
    model = lm( y ~ grs + I( grs^2 ))
    init = model$coefficients[ 2:3 ]
    var = sum((model$residuals)^2)/(n-3)
    coeff = .single_gxe(y,grs,c(init, 0.1, 0, var))
  } else if (method=='OLS'){
    x1 = cbind(grs,(grs^2-1))
    par1 = OLS_1(y,x1,b0=F)
    x2 = cbind(grs,grs^2)
    # yr2 = (y - x1%*%par1)^2
    # copar = lm(yr2~x2)$coefficients
    # par2 = nls(yr2 ~ (b+c*grs)^2+sig^2,start=list(b=0.1,c=0,sig=0.5))
    # par2 = OLS_1(yr2,x2,b0=T)
    # gamma = sqrt(abs(par2[3]))
    # beta = par2[2]/(2*gamma)
    # sig = sqrt(par2[1] - beta^2)
    # coeff = c(par1,beta,gamma,sig)
    coeff = OLS_2(y,x2,par1)
  }
  coef_names = c('alpha1', 'alpha2', 'beta', 'gamma','sigma' )
  return(setNames( coeff, coef_names )) 
}

################ simulation ###################
m    = 100       # number of genetic markers
n    = 1e4       # sample size
a0   = sqrt(.1)  # linear effect of GRS on y
b1   = sqrt(.3)  # linear effect of E on y
c1   = sqrt(.05)  # interaction effect
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E
skwN = 0         # skewness of the noise
krtN = 3         # kurtosis of the noise
pow  = 1         # transformation power

d   = seq(0,0.3,0.05)         # correlation between E and GRS
# d=0.3

# m    = 100       # number of genetic markers
# n    = 1e4       # sample size
# a1   = sqrt(.1)  # linear effect of GRS on y
# a2   = 0
# b1   = sqrt(.3)  # beta
# # sigE = sqrt(.5)
# # c1   = sqrt(.05)  # gamma
# sig  = 0.5
# skwE = 0         # skewness of E
# krtE = 3         # kurtosis of E
# skwN = 0         # skewness of the noise
# krtN = 3         # kurtosis of the noise
# pow  = 1
# 
# c = sqrt(seq(0,0.3,0.05))

if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
  stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
}

# Simulation

cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)
# clusterExport(cl,c(),envir = environment())
Sys.time()

basic_twostep_OLS_m1 = foreach(d1=d,.export=c('rpearson')) %dopar% {
  # start=Sys.time()
  # alld = lapply(t(as.matrix(d)),function(d1){
  reps = sapply(1:500,function(rep){
    maf  =  matrix( runif( m ), nrow = 1 )
    G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    eff  =  matrix( runif( m ), ncol = 1 )
    eff  =  eff / sqrt( sum( eff^2 ) )
    GRS  =  scale( G %*% eff )

    E    =  d1 * GRS + sqrt( 1 - d1^2 ) * matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
    noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
    sig  =  sqrt( 1 - a0^2 - b1^2 - c1^2 * (1 + d1^2) - 2 * a0 * b1 * d1 ) ###transformation of sig

    z    = a0 * GRS + b1 * E + c1 * GRS * E + sig * noi
    y    =  z - mean(z)
    
    # maf  =  matrix( runif( m ), nrow = 1 )
    # G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    # eff  =  matrix( runif( m ), ncol = 1 )
    # eff  =  eff / sqrt( sum( eff^2 ) )
    # GRS  =  scale( G %*% eff )
    # 
    # E    =  matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
    # noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
    # # sig  =  sqrt( 1 - a1^2 - 2*a2^2 - b1^2 - c1^2 )   ###transformation of sig
    # 
    # # z    =  scale( a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi )
    # z    =  a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi 
    # y = z - mean(z)
    
    # Ftrans  =  function( s, p1, p2 ) ( (s-p1)^p2 - 1 ) / p2 
    # 
    # if (pow != 0) {
    #   # y  =  scale( Ftrans( z, min(z)-1e-5, pow ) )
    #   y = Ftrans( z, min(z)-1e-5, pow )
    # } else {
    #   y  =  scale( log( z - min(z)+1e-5 ) )
    # }
    
    # sel  =  which( abs( y ) > 10 )
    # 
    # while (length( sel ) > 0) {
    #   y[ sel ] = NA
    #   y  =  scale( y )
    #   sel  =  which( abs( y ) > 10 )
    # }
    # 
    # idx = which(is.na(y))
    # if (length(idx) != 0){
    #   y = y[-idx]
    #   GRS= as.matrix(GRS[-idx,])
    # } else {
    #   y = y
    #   GRS = GRS
    # }
    
    # Estimate interaction effect for GRS
    return( estimate_basic( y, GRS, method='OLS' ) )
  })
  # Sys.time()-start
  return(reps)
}
stopCluster(cl)
Sys.time()
save(basic_twostep_OLS_m1,file='./data/new/basic_twostep_OLS_m1.Rdata')

data = basic_twostep_OLS_m1
# alpha0 = sapply(basic_twostep,function(x) x['alpha0',])
alpha1 = sapply(data,function(x) x['alpha1',])
alpha2 = sapply(data,function(x) x['alpha2',])
beta = sapply(data,function(x) x['beta',])
gamma = sapply(data,function(x) x['gamma',])
sigma = sapply(data,function(x) x['sigma',])
# sigmaE = sapply(basic,function(x) x['sigmaE',])
colnames(alpha1) <- colnames(alpha2) <- colnames(beta) <- colnames(gamma) <- colnames(sigma) <- d

gamma0 = sqrt(alpha2^2+gamma^2)
delta = alpha2/gamma0
beta0 = beta/sqrt(1-delta^2)
alpha0 = alpha1-beta0*delta

delta[delta>0.9]=NA

# pdf('./results/new/basic_twostep_OLS.pdf',width = 10,height=12)
# par(mfrow=c(3,2))
# # boxplot(alpha0,xlab='different gamma^2',ylab='Estimation for a0',las=1)
# # abline(h=a0,lty=2)
# boxplot(alpha1,xlab='different gamma^2',ylab='Estimation for alpha1',main='alpha1',las=1)
# abline(h=a1,lty=2)
# boxplot(alpha2,xlab='different gamma^2',ylab='Estimation for alpha2',main='alpha2',las=1)
# abline(h=a2,lty=2)
# boxplot(beta,xlab='different gamma^2',ylab='Estimation for beta',main='beta',las=1)
# abline(h=b1,lty=2)
# boxplot(gamma^2,xlab='different gamma^2',ylab='Estimation for gamma^2',main='gamma',las=1)
# arrows(c(1:7)-0.6,c^2,c(1:7)+0.6,c^2,angle=0,lty=2)
# boxplot(sigma,xlab='different gamma^2',ylab='Estimation for sigma',main='sigma',las=1)
# abline(h=sig,lty=2)
# dev.off()

pdf('./results/new/basic_twostep_OLS_m1.pdf',width = 10,height=8)
par(mfrow=c(2,2))
boxplot(alpha0,xlab='different delta',ylab='Estimation for alpha0',main='alpha1',las=1)
abline(h=a0,lty=2)
boxplot(beta0,xlab='different delta',ylab='Estimation for beta0',main='beta',las=1)
abline(h=b1,lty=2)
boxplot(gamma0,xlab='different delta',ylab='Estimation for gamma0',main='gamma',las=1)
abline(h=c1,lty=2)
boxplot(delta,xlab='different delta',ylab='Estimation for delta',main='delta',las=1)
arrows(c(1:7)-0.6,d,c(1:7)+0.6,d,angle=0,lty=2)
dev.off()

# pdf('./results/new/basic_twostep_OLS_gamma.pdf',width = 10,height=4)
# par(mfrow=c(1,2))
# boxplot(beta[,-1],xlab='different gamma^2',ylab='Estimation for beta',main='supplementary detail of beta',las=1)
# abline(h=b1,lty=2)
# boxplot(sigma[,-1],xlab='different gamma^2',ylab='Estimation for sigma',main='supplementary detail of sigma',las=1)
# abline(h=sig,lty=2)
# dev.off()

# pdf('./results/new/basic.pdf',width = 10,height=8)
# par(mfrow=c(2,2))
# boxplot(alpha0,xlab='delta',ylab='Estimation for alpha',las=1)
# abline(h=a1,lty=2)
# boxplot(beta0,xlab='delta',ylab='Estimation for beta',las=1)
# abline(h=b1,lty=2)
# boxplot(gamma0,xlab='delta',ylab='Estimation for gamma',las=1)
# abline(h=c1,lty=2)
# boxplot(delta,xlab='delta',ylab='Estimation for delta',yaxt='n')
# axis(2,at=seq(-0.1,0.4,0.05),las=1)
# arrows(c(1:7)-0.6,d,c(1:7)+0.6,d,angle=0,lty=2)
# dev.off()



################ skw krt change ##################
m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
c1   = sqrt(0.025)  # sqrt(.3)  # gamma
sig  = 0.5
# skwE = 0         # skewness of E
# krtE = 3         # kurtosis of E
# allskwN = c(0,2,4) # skewness of the noise
# allkrtN = c(2,6,8,18) # kurtosis of the noise
pow  = 1         # transformation power

# skwN = 0         # skewness of E
# krtN = 3         # kurtosis of E
# allskwE = c(0,1,2,3,4,5) # skewness of the noise
# krtE = 27
# # skwE = 0
# # allkrtE = c(2,3,6,11,18,27) # kurtosis of the noise

skwE = 0         # skewness of E
krtE = 3
# allskwN = c(0,1,2,3,4,5) # skewness of the noise
# krtN = 27
skwN = 0
allkrtN = c(2,3,6,11,18,27) # kurtosis of the noise

# allcomb = expand.grid(allskwE,allkrtE)
# allcomb = allcomb[which(allcomb[,2] > allcomb[,1]^2+1),]
# colnames(allcomb) = c('skwE','krtE')
# allcomb = t(allcomb)

# Simulation

cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)
# clusterExport(cl,c(),envir = environment())

Sys.time()
# basic_twostep_OLS_skwE <- foreach(comb=1:ncol(allcomb),.export=c('rpearson')) %dopar% {
basic_twostep_OLS_ind_krtN <- foreach(krtN=allkrtN,.export=c('rpearson')) %dopar% {
    # sim21 = parLapply(cl,allcomb,function(x){
  # skwE = allcomb[1,comb]
  # krtE = allcomb[2,comb]
  
  reps = sapply(1:500,function(rep){
    if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
      stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
    }
    
    maf  =  matrix( runif( m ), nrow = 1 )
    G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    eff  =  matrix( runif( m ), ncol = 1 )
    eff  =  eff / sqrt( sum( eff^2 ) )
    GRS  =  scale( G %*% eff )
    
    E    =  matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
    noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
    # sig  =  sqrt( 1 - a1^2 - 2*a2^2 - b1^2 - c1^2 )   ###transformation of sig
    
    # z    =  scale( a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi )
    z    =  a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi 
    
    y = z - mean(z)
    
    # Estimate interaction effect for GRS
    return( estimate_basic( y, GRS,method='OLS' ) )
  })
  # Sys.time()-start
  return(reps)
}
stopCluster(cl)
Sys.time()
save(basic_twostep_OLS_ind_krtN,file='./data/new/basic_twostep_OLS_ind_krtN.Rdata')

# data = basic_twostep_OLS_ind_krtE
# # a0 = sapply(basic_twostep_skwE,function(x) x['a0',])
# alpha1 = sapply(data,function(x) x['alpha1',])
# alpha2 = sapply(data,function(x) x['alpha2',])
# beta = sapply(data,function(x) x['beta',])
# gamma = sapply(data,function(x) x['gamma',])
# sigma = sapply(data,function(x) x['sigma',])
# colnames(alpha1) <- colnames(alpha2) <- colnames(beta) <- colnames(gamma) <- colnames(sigma) <- allkrtE

# beta0 = b1^2 + sig^2
# beta1 = 2*b1*c1
# beta2 = c1^2

# pdf('./results/new/basic_twostep_OLS_ind_skwE.pdf',width = 10,height=8)
# par(mfrow=c(3,2))
# boxplot(alpha1,xlab='skw and krt',ylab='Estimation for alpha1',main='alpha1',las=1)
# abline(h=a1,lty=2)
# boxplot(alpha2,xlab='skw and krt',ylab='Estimation for alpha2',,main='alpha2',las=1)
# abline(h=a2,lty=2)
# boxplot(beta,xlab='skw and krt',ylab='Estimation for beta',main='beta',las=1)
# abline(h=b1,lty=2)
# boxplot(gamma,xlab='skw and krt',ylab='Estimation for gamma',main='gamma',las=1)
# abline(h=c1,lty=2)
# boxplot(sigma2,xlab='skw and krt',ylab='Estimation for sigma^2',main='sigma^2',las=1)
# abline(h=sig^2,lty=2)
# dev.off()

skwE_g = sapply(basic_twostep_OLS_ind_skwE,function(x) x['gamma',])
krtE_g = sapply(basic_twostep_OLS_ind_krtE,function(x) x['gamma',])

s.var=apply(skwE_g,2,var)
s.bias = apply(skwE_g,2,mean)-c1
k.var=apply(krtE_g,2,var)
k.bias = apply(krtE_g,2,mean)-c1

pdf('./results/new/basic_twostep_OLS_indE_bias_var.pdf',width = 8,height=8)
par(mfrow=c(2,2))
plot(s.bias~allskwE,type='l',xlab='skweness(E)',las=1)
plot(s.var~allskwE,type='l',xlab='skweness(E)',ylab='Variance',las=1)
plot(k.bias~allkrtE,type='l',xlab='kurtosis(E)',las=1,xaxt='n')
axis(side=1,at=allkrtE)
plot(k.var~allkrtE,type='l',xlab='kurtosis(E)',ylab='Variance',las=1,xaxt='n')
axis(side=1,at=allkrtE)
dev.off()


skwN_g = sapply(basic_twostep_OLS_ind_skwN,function(x) x['gamma',])
krtN_g = sapply(basic_twostep_OLS_ind_krtN,function(x) x['gamma',])

s.var=apply(skwN_g,2,var)
s.bias = apply(skwN_g,2,mean)-c1
k.var=apply(krtN_g,2,var)
k.bias = apply(krtN_g,2,mean)-c1

pdf('./results/new/basic_twostep_OLS_indN_bias_var.pdf',width = 8,height=8)
par(mfrow=c(2,2))
plot(s.bias~allskwN,type='l',xlab='skweness(N)',las=1)
plot(s.var~allskwN,type='l',xlab='skweness(N)',ylab='Variance',las=1)
plot(k.bias~allkrtN,type='l',xlab='kurtosis(N)',las=1,xaxt='n')
axis(side=1,at=allkrtN)
plot(k.var~allkrtN,type='l',xlab='kurtosis(N)',ylab='Variance',las=1,xaxt='n')
axis(side=1,at=allkrtN)
dev.off()

# parname = c('alpha1','alpha2','beta','gamma','sigma')
# parvalue = c(a1,a2,b1,c1,sig)
# names(parvalue) = parname
# 
# plot=list()
# n = 0
# for (par in parname){
#   data = sapply(basic_twostep_OLS_skwE,function(x) x[par,])
#   data = cbind(t(allcomb),t(data))
#   data = as.data.frame(data)
#   data[,1] = factor(data[,1])
#   data[,2] = factor(data[,2],allkrtE)
#   colnames(data) = c('skwE','krtE',1:500)
#   datas = melt(data[,-2],id.vars=c('skwE'),varnames ='rep',value.name = 'data')
#   datak = melt(data[,-1],id.vars=c('krtE'),varnames ='rep',value.name = 'data')
#   
#   # if (par=='gamma'){
#   #   datas = datas[-which(datas$data>0.9),]
#   #   datak = datak[-which(datak$data>0.9),] 
#   # }
#   
#   plot[[n+1]]= ggplot(datas,aes(skwE,data))+geom_boxplot()+
#     geom_hline(yintercept = parvalue[par],linetype='dashed',color='red')+
#     labs(x='Skewness(E)',y=par,title='gamma^2=0.025')
#     # lims(y=c(0.4,0.7))
#   plot[[n+2]]= ggplot(datak,aes(krtE,data))+geom_boxplot()+
#     geom_hline(yintercept = parvalue[par],linetype='dashed',color='red')+
#     labs(x='Kurtosis(E)',y=par ,title='gamma^2=0.025')
#     # lims(y=c(0.1,0.2))
#   n = n+2
# }

# # par='gamma'
# # data = sapply(basic_twostep_OLS_skwE_g0,function(x) x[par,])
# # data2 = melt(data,id.vars=c('skwE','krtE'),varnames ='rep',value.name = 'value')
# # datass = dcast(data2,krtE+variable~skwE,value.var='value')
# # datakk = dcast(data2,skwE+variable~krtE,value.var='value')
# # 
# # pdf('./results/new/basic_twostep_OLS_skwE_g0_gamma.pdf',width = 10,height=5)
# # par(mfrow=c(1,2))
# # boxplot(datass[,-c(1:2)],outline = F,las=1,xlab='Skewness(E)',ylab=par,main='gamma=0')
# # boxplot(datakk[,-c(1:2)],outline = F,las=1,xlab='Kurtosis(E)',ylab=par,main='gamma=0')
# # dev.off()
# 

pdf('./results/new/basic_twostep_OLS_skwE.pdf',width = 8,height=15)
grid.arrange(grobs=plot,nrow=5)
dev.off()

pdf('./results/new/basic_twostep_OLS_skwE_gamma.pdf',width = 8,height=3)
grid.arrange(grobs=c(plot[7:8]),nrow=1)
dev.off()

pdf('./results/new/basic_twostep_OLS_skwE_beta.pdf',width = 8,height=6)
grid.arrange(grobs=c(plot[5:6],plot[11:12]),nrow=2)
dev.off()

# out = data2[which(data2$value>0.9),]




################ N skw krt change ##################
m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
c1   = sqrt(.025)  # gamma
sig  = 0.5
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E
allskwN = c(0,1,2,3,4,5) # skewness of the noise
allkrtN = c(2,3,6,11,18,27) # kurtosis of the noise
pow  = 1         # transformation power

# skwN = 0         # skewness of E
# krtN = 3         # kurtosis of E
# allskwE = c(0,1,2,3,4,5) # skewness of the noise
# allkrtE = c(2,3,6,11,18,27) # kurtosis of the noise

allcomb = expand.grid(allskwN,allkrtN)
allcomb = allcomb[which(allcomb[,2] > allcomb[,1]^2+1),]
colnames(allcomb) = c('skwN','krtN')
allcomb = t(allcomb)

# Simulation

cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)
# clusterExport(cl,c(),envir = environment())

Sys.time()
basic_twostep_OLS_skwN <- foreach(comb=1:ncol(allcomb),.export=c('rpearson')) %dopar% {
  # sim21 = parLapply(cl,allcomb,function(x){
  skwN = allcomb[1,comb]
  krtN = allcomb[2,comb]
  
  reps = sapply(1:500,function(rep){
    if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
      stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
    }
    
    maf  =  matrix( runif( m ), nrow = 1 )
    G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    eff  =  matrix( runif( m ), ncol = 1 )
    eff  =  eff / sqrt( sum( eff^2 ) )
    GRS  =  scale( G %*% eff )
    
    E    =  matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
    noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
    # sig  =  sqrt( 1 - a1^2 - 2*a2^2 - b1^2 - c1^2 )   ###transformation of sig
    
    # z    =  scale( a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi )
    z    =  a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi 
    
    y = z - mean(z)
    
    # Estimate interaction effect for GRS
    return( estimate_basic( y, GRS, method='OLS' ) )
  })
  # Sys.time()-start
  return(reps)
}
stopCluster(cl)
Sys.time()
save(basic_twostep_OLS_skwN,file='./data/new/basic_twostep_OLS_skwN.Rdata')

# data = basic_twostep_OLS_skwN_g0
# # a0 = sapply(basic_twostep_skwE,function(x) x['a0',])
# alpha1 = sapply(data,function(x) x['alpha1',])
# alpha2 = sapply(data,function(x) x['alpha2',])
# beta = sapply(data,function(x) x['beta',])
# gamma = sapply(data,function(x) x['gamma',])
# sigma = sapply(data,function(x) x['sigma',])

# pdf('./results/new/basic_twostep_skwE.pdf',width = 10,height=8)
# par(mfrow=c(2,2))
# boxplot(alpha1,xlab='skw and krt',ylab='Estimation for alpha1',main='alpha1',las=1)
# abline(h=a1,lty=2)
# boxplot(alpha2,xlab='skw and krt',ylab='Estimation for alpha2',,main='alpha2',las=1)
# abline(h=a2,lty=2)
# boxplot(beta,xlab='skw and krt',ylab='Estimation for beta',main='beta',las=1)
# abline(h=b1,lty=2)
# boxplot(gamma,xlab='skw and krt',ylab='Estimation for gamma',main='gamma',las=1)
# abline(h=c1,lty=2)
# boxplot(sigma,xlab='skw and krt',ylab='Estimation for sigma',main='sigma',las=1)
# abline(h=sig,lty=2)
# dev.off()

parname = c('alpha1','alpha2','beta','gamma','sigma')
parvalue = c(a1,a2,b1,c1,sig)
names(parvalue) = parname

plot=list()
n = 0
for (par in parname){
  data = sapply(basic_twostep_OLS_skwN,function(x) x[par,])
  data = cbind(t(allcomb),t(data))
  data = as.data.frame(data)
  data[,1] = factor(data[,1])
  data[,2] = factor(data[,2],allkrtN)
  colnames(data) = c('skwN','krtN',1:100)
  datas = melt(data[,-2],id.vars=c('skwN'),varnames ='rep',value.name = 'data')
  datak = melt(data[,-1],id.vars=c('krtN'),varnames ='rep',value.name = 'data')
  
  plot[[n+1]]= ggplot(datas,aes(skwN,data))+geom_boxplot()+
    geom_hline(yintercept = parvalue[par],linetype='dashed',color='red')+
    labs(x='Skewness(N)',y=par,title='gamma^2=0.025')+
    lims(y=c(0.4,0.7))
  plot[[n+2]]= ggplot(datak,aes(krtN,data))+geom_boxplot()+
    geom_hline(yintercept = parvalue[par],linetype='dashed',color='red')+
    labs(x='Kurtosis(N)',y=par,title='gamma^2=0.025')+
    lims(y=c(0.4,0.7))
  n = n+2
}


pdf('./results/new/basic_twostep_OLS_skwN.pdf',width = 8,height=15)
grid.arrange(grobs=plot,nrow=5)
dev.off()

pdf('./results/new/basic_twostep_OLS_skwN_gamma.pdf',width = 8,height=3)
grid.arrange(grobs=c(plot[7:8]),nrow=1)
dev.off()

pdf('./results/new/basic_twostep_OLS_skwN_beta.pdf',width = 8,height=6)
grid.arrange(grobs=c(plot[5:6],plot[11:12]),nrow=2)
dev.off()


#### test for var(b)
m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
c1   = sqrt(0.025)  # sqrt(.3)  # gamma
sig  = 0.5
# skwE = 0         # skewness of E
# krtE = 3         # kurtosis of E
# allskwN = c(0,2,4) # skewness of the noise
# allkrtN = c(2,6,8,18) # kurtosis of the noise
pow  = 1         # transformation power

skwN = 0         # skewness of E
krtN = 3         # kurtosis of E
allskwE = c(0,1,2,3,4,5) # skewness of the noise
allkrtE = c(2,3,6,11,18,27) # kurtosis of the noise

allcomb = expand.grid(allskwE,allkrtE)
allcomb = allcomb[which(allcomb[,2] > allcomb[,1]^2+1),]
colnames(allcomb) = c('skwE','krtE')
allcomb = t(allcomb)


# Simulation

cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)
# clusterExport(cl,c(),envir = environment())

Sys.time()
test_var <- foreach(comb=1:ncol(allcomb),.export=c('rpearson')) %dopar% {
  # sim21 = parLapply(cl,allcomb,function(x){
  skwE = allcomb[1,comb]
  krtE = allcomb[2,comb]
  
  if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
    stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
  }
  
  maf  =  matrix( runif( m ), nrow = 1 )
  G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
  eff  =  matrix( runif( m ), ncol = 1 )
  eff  =  eff / sqrt( sum( eff^2 ) )
  GRS  =  scale( G %*% eff )
  
  E    =  matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
  noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
  # sig  =  sqrt( 1 - a1^2 - 2*a2^2 - b1^2 - c1^2 )   ###transformation of sig
  
  # z    =  scale( a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi )
  z    =  a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi 
  
  y = z - mean(z)
  
  # calculate var for GRS
  grs = GRS
  x1 = cbind(grs,(grs^2-1))
  par1 = OLS_1(y,x1,b0=F)
  x2 = cbind(grs,grs^2)
  y2 = (y - x1%*%par1)^2
  par2 = OLS_1(y2,x2,b0=T)
  var.s1 = var(y - x1%*%par1)
  n = nrow(x2)
  p = ncol(x2)
  sig22 = sum((y2-x2%*%par2[2:3]-par2[1])^2)/(n-p-1)
  cov = calcov(y2,x2)
  
  # var.gamma2 = cov[3,3]
  return(list(cov=cov,var.s1=var.s1,mean.y2=mean(y2),sig22=sig22))
}
stopCluster(cl)

data = sapply(test_var,function(x) x[[4]])
data = cbind(t(allcomb),data=data)
# data = mean.y2
aggregate(data~skwE,data=data,mean)
aggregate(data~krtE,data=data,mean)
boxplot(mean~skwE,data=data)
boxplot(mean~krtE,data=data)
## change in sig2 !!!

calcov = function(y2,x2){
  n = nrow(x2)
  p = ncol(x2)
  
  x2 = cbind(1,x2)
  
  # xba = colMeans(x2)
  # one = rep(1,n)
  # xres = x2 - one%*%t(xba)
  # 
  # xsol = solve(t(xres)%*%xres)
  # c = one %*% t(xba) %*% xsol %*% t(xres) 
  # cov = sum((c %*% y2)^2) / (n-p-1) * xsol
  
  xsol = solve(t(x2)%*%x2)
  b = xsol %*% t(x2) %*% y2
  cov = sum((y2-x2%*%b)^2) / (n-p-1) * xsol
  
  cov
}
