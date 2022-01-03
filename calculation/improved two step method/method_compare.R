setwd('D:/Jiang Lab/GxE')

library(sjPlot)

### normal distribution
load('./data/new/basic_twostep_OLS_m1.Rdata') ## data: basic_twostep_OLS_m1
load('./data/fig1.Rdata')   ##data: sim1

## no skw krt
a0   = sqrt(.1)  # linear effect of GRS on y
b1   = sqrt(.3)  # linear effect of E on y
c1   = sqrt(.05)  # interaction effect
d   = seq(0,0.3,0.05)         # correlation between E and GRS

alpha1 = a0 + b1*d
alpha2 = c1*d
beta   = b1*sqrt(1-d^2)
gamma  = c1*sqrt(1-d^2)
sig    = sqrt( 1 - a0^2 - b1^2 - c1^2 * (1 + d^2) - 2 * a0 * b1 * d )

data = basic_twostep_OLS_m1
parname = c('alpha1','alpha2','beta','gamma') # ,'sigma')
parvalue = cbind(alpha1,alpha2,beta,gamma)  #,sig)
colnames(parvalue) = parname

normal = NULL
for (par in parname){
  par_new = sapply(data,function(x) x[par,])
  bias_n = colMeans(par_new)-parvalue[,par]
  var_n  = apply(par_new,2,var)
  # sd_n   = apply(par_new,2,sd)
  # try1 = rbind(bias_n,var_n)
  # try1 = as.numeric(try1)
  # mean_bias = mean(bias_n)
  try1 = paste0(round(bias_n,4),'(',round(var_n,4),')')
  
  par_old = lapply(sim1,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][[par]]))
  par_old = do.call(rbind,par_old)
  bias_o = colMeans(par_old)-parvalue[,par]
  var_o  = apply(par_old,2,var)
  # try2 = rbind(bias_o,var_o)
  # try2 = as.numeric(try2)
  try2 = paste0(round(bias_o,4),'(',round(var_o,4),')')
  
  try = rbind(try1,try2)
  try = cbind(par=par,method=c('OLS','Article'),as.data.frame(try))
  normal = rbind(normal,try)
}

# idx1 = rbind(paste0(d,'.bias'),paste0(d,'.var'))
# dim(idx1) = dim(idx1)[1]*dim(idx1)[2]
colnames(normal) = c('par','method',d)

tab_df(normal)



### skwE krtE
load('./data/fig2.Rdata')  ## sim21: skwN   sim22: skwE
load('./data/new/basic_twostep_OLS_skwE.Rdata')

a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
c1   = sqrt(.025)         # gamma
allskwE = c(0,2,4) # skewness of the noise
allkrtE = c(2,6,18) # kurtosis of the noise
allcomb = expand.grid(allskwE,allkrtE)
allcomb = allcomb[which(allcomb[,2] > allcomb[,1]^2+1),]
colnames(allcomb) = c('skwE','krtE')
allcomb = t(allcomb)

# allskwE = c(0,1,2,3,4,5) # skewness of the noise
# allkrtE = c(2,3,6,11,18,27) # kurtosis of the noise
# the selected columns: 1,4,6,11,13,15
select_OLS_E = basic_twostep_OLS_skwE[c(1,4,6,11,13,15)]

parname = c('alpha1','alpha2','beta','gamma') # ,'sigma')
parvalue = cbind(a1,a2,b1,c1)  #,sig)
parvalue2 = cbind(a1,a2,b1,0)
names(parvalue) <- names(parvalue2) <- parname

E = NULL
for (par in parname){
  data = sapply(select_OLS_E,function(x) x[par,])
  data = cbind(t(allcomb),t(data[1:100,]))
  data = as.data.frame(data)
  data[,1] = factor(data[,1])
  data[,2] = factor(data[,2],allkrtE)
  colnames(data) = c('skwE','krtE',1:100)
  datas = melt(data[,-2],id.vars=c('skwE'),varnames ='rep',value.name = 'data')
  datak = melt(data[,-1],id.vars=c('krtE'),varnames ='rep',value.name = 'data')
  bias_ns = aggregate(data~skwE,data=datas,mean)[,2] - parvalue[par]
  var_ns = aggregate(data~skwE,data=datas,'var')[,2]
  bias_nk = aggregate(data~krtE,data=datak,mean)[,2] - parvalue[par]
  var_nk = aggregate(data~krtE,data=datak,'var')[,2]
  
  par_old = lapply(sim22,function(x) c(x$par.skwN,x$par.krtN, sapply(x[1:100], function(n) n[["real_data"]][["coefficients"]][[par]])))
  par_old = do.call(rbind,par_old)
  par_old = as.data.frame(par_old)
  par_old[,1] = factor(par_old[,1])
  par_old[,2] = factor(par_old[,2],allkrtE)
  colnames(par_old) = c('skwE','krtE',1:100)
  par_s = melt(par_old[,-2],id.vars=c('skwE'),varnames ='rep',value.name = 'data')
  par_k = melt(par_old[,-1],id.vars=c('krtE'),varnames ='rep',value.name = 'data')
  bias_os = aggregate(data~skwE,data=par_s,mean)[,2] - parvalue2[par]
  var_os = aggregate(data~skwE,data=par_s,'var')[,2]
  bias_ok = aggregate(data~krtE,data=par_k,mean)[,2] - parvalue2[par]
  var_ok = aggregate(data~krtE,data=par_k,'var')[,2]
  
  bias_n = round(c(bias_ns,bias_nk),5)
  var_n  = round(c(var_ns,var_nk),5)
  new = paste0(bias_n,'(',var_n,')')
    # as.numeric(rbind(c(bias_ns,bias_nk),c(var_ns,var_nk)))
  bias_o = round(c(bias_os,bias_ok),5)
  var_o  = round(c(var_os,var_ok),5)
  old = paste0(bias_o,'(',var_o,')')
  try = rbind(new,old)
  try = cbind(par=par,method=c('OLS','Article'),as.data.frame(try))
  E = rbind(E,try)
}

# idx1 = rep(c('bias','var'),6)
# idx2 = c(rep(paste0('skwE=',allskwE),each=2),rep(paste0('krtE_',allkrtE),each=2))
colnames(E) = c('par','method',paste0('skwE=',allskwE),paste0('krtE=',allkrtE))
# E[,3:14] = round(E[,3:14],5)




### skwN krtN
load('./data/fig2.Rdata')  ## sim21: skwN   sim22: skwE
load('./data/new/basic_twostep_OLS_skwN.Rdata')

a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
c1   = sqrt(.025)        # gamma
allskwN = c(0,2,4) # skewness of the noise
allkrtN = c(2,6,18) # kurtosis of the noise
allcomb = expand.grid(allskwN,allkrtN)
allcomb = allcomb[which(allcomb[,2] > allcomb[,1]^2+1),]
colnames(allcomb) = c('skwN','krtN')
allcomb = t(allcomb)

# allskwN = c(0,1,2,3,4,5) # skewness of the noise
# allkrtN = c(2,3,6,11,18,27) # kurtosis of the noise
# the selected columns: 1,4,6,11,13,15
select_OLS_N = basic_twostep_OLS_skwN[c(1,4,6,11,13,15)]

parname = c('alpha1','alpha2','beta','gamma') # ,'sigma')
parvalue = cbind(a1,a2,b1,c1)  #,sig)
parvalue2 = cbind(a1,a2,b1,0)
names(parvalue) <- names(parvalue2) <- parname

N = NULL
for (par in parname){
  data = sapply(select_OLS_N,function(x) x[par,])
  data = cbind(t(allcomb),t(data[1:100,]))
  data = as.data.frame(data)
  data[,1] = factor(data[,1])
  data[,2] = factor(data[,2],allkrtN)
  colnames(data) = c('skwN','krtN',1:100)
  datas = melt(data[,-2],id.vars=c('skwN'),varnames ='rep',value.name = 'data')
  datak = melt(data[,-1],id.vars=c('krtN'),varnames ='rep',value.name = 'data')
  bias_ns = aggregate(data~skwN,data=datas,mean)[,2] - parvalue[par]
  var_ns = aggregate(data~skwN,data=datas,'var')[,2]
  bias_nk = aggregate(data~krtN,data=datak,mean)[,2] - parvalue[par]
  var_nk = aggregate(data~krtN,data=datak,'var')[,2]
  
  par_old = lapply(sim21,function(x) sapply(x[1:100], function(n) n[["real_data"]][["coefficients"]][[par]]))
  par_old = do.call(rbind,par_old)
  par_old = cbind(t(allcomb),par_old)
  par_old = as.data.frame(par_old)
  par_old[,1] = factor(par_old[,1])
  par_old[,2] = factor(par_old[,2],allkrtN)
  colnames(par_old) = c('skwN','krtN',1:100)
  par_s = melt(par_old[,-2],id.vars=c('skwN'),varnames ='rep',value.name = 'data')
  par_k = melt(par_old[,-1],id.vars=c('krtN'),varnames ='rep',value.name = 'data')
  bias_os = aggregate(data~skwN,data=par_s,mean)[,2] - parvalue2[par]
  var_os = aggregate(data~skwN,data=par_s,'var')[,2]
  bias_ok = aggregate(data~krtN,data=par_k,mean)[,2] - parvalue2[par]
  var_ok = aggregate(data~krtN,data=par_k,'var')[,2]
  
  bias_n = round(c(bias_ns,bias_nk),5)
  var_n  = round(c(var_ns,var_nk),5)
  new = paste0(bias_n,'(',var_n,')')
  # as.numeric(rbind(c(bias_ns,bias_nk),c(var_ns,var_nk)))
  bias_o = round(c(bias_os,bias_ok),5)
  var_o  = round(c(var_os,var_ok),5)
  old = paste0(bias_o,'(',var_o,')')
  try = rbind(new,old)
  try = cbind(par=par,method=c('OLS','Article'),as.data.frame(try))
  N = rbind(N,try)
}

# idx1 = rep(c('bias','var'),6)
# idx2 = c(rep(paste0('skwN_',allskwN),each=2),rep(paste0('krtN_',allkrtN),each=2))
colnames(N) = c('par','method',paste0('skwN=',allskwN),paste0('krtN=',allkrtN))
# N[,3:14] = round(N[,3:14],5)

tab_df(E, title='Skwness / Kurtosis of E')
tab_df(N, title='Skwness / Kurtosis of Noise')

### for gamma
Eg = E[7:8,]
Ng = N[7:8,]
colnames(Eg) <- colnames(Ng) <-  c('par','method',paste0('skw=',allskwN),paste0('krt=',allkrtN))
gamma = t(cbind(variable=rep(c('E','Noise'),each=2),rbind(Eg,Ng)))
gamma = as.data.frame(gamma)
gamma = cbind(rownames(gamma),gamma)
# colnames(gamma) = paste(gamma[1,],gamma[2,],gamma[3,],sep='_')
tab_df(gamma,title='Estimation of gamma')
