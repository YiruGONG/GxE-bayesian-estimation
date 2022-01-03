computer = 1    ### 1 for own computer ; 0 for thinkpad
dirct = ifelse(computer,'D:/Jiang Lab/','E:/')
setwd(paste0(dirct,'GxE'))

# library('devtools')
# install_github('zkutalik/GRSxE_software',subdir='Rcode')
# library( GxE )
library( PearsonDS )
library(RNOmni)
library(doParallel)
library(foreach)
library(gsl)

source('./GRSxE_software/Rcode/R/estimate_gxe.R')


############ supp - Outcome modelled on transformed scale ############
## n = 10000, a1^2=0.1,a2=0, beta^2=0.3, E~N(0,1),f(t)=t
##  gamma^2=0 or gamma^2=0.05
## transformation f(t)=t^p, p=-2,-1
## generate fake GRS

m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.05)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E
# pow  = 1
# c1 = sqrt(.05)
# skwN = 0
# krtN = 2

allskwN = c(0,1,2)    # skewness of the noise
allkrtN = c(2,3,6,8) # kurtosis of the noise
# allc1   = sqrt(0.05)        # gamma
allc1 = c(0,sqrt(0.05))
allpow  = 11  #'xlogx' # transformation power

# allcomb = expand.grid(allc1,allskwN,allkrtN)
# allcomb = allcomb[which(allcomb[,3] > allcomb[,2]^2+1),]
# colnames(allcomb) = c('c1','skwN','krtN')
# allcomb = t(allcomb)

allcomb = expand.grid(allc1,allpow,allskwN,allkrtN)
# allcomb = allcomb[-c(1,2,13,14,17,18),]
allcomb = allcomb[which(allcomb[,4] > allcomb[,3]^2+1),]
colnames(allcomb) = c('c1','power','skwN','krtN')
allcomb = t(allcomb)

# Simulation
###open multiple cores
cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)

sim_xlogx_nor <- foreach(rep=1:50,.export=c('rpearson','estimate_gxe','mclapply')) %dopar% {
  all <- lapply(as.data.frame(allcomb),function(x){
    c1 = x[1]
    pow = x[2]
    skwN = x[3]
    krtN = x[4]

    if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
      stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
    }
    
    sim_num = 1e2 # number of bootstrap / fake GRS
    
    maf  =  matrix( runif( m ), nrow = 1 )
    G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    eff  =  matrix( runif( m ), ncol = 1 )
    eff  =  eff / sqrt( sum( eff^2 ) )
    GRS  =  scale( G %*% eff )
    
    E    =  matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
    noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
    sig  =  sqrt( 1 - a1^2 - 2*a2^2 - b1^2 - c1^2 )   ###transformation of sig
    
    z    =  scale( a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi )
    
    # Ftrans  =  function( s, p1, p2 ) ( (s-p1)^p2 - 1 ) / p2
    # 
    # if (pow != 0) {
    #   y  =  scale( Ftrans( z, min(z)-1e-5, pow ) )
    # } else {
    #   y  =  scale( log( z - min(z)+1e-5 ) )
    # }
    Ftrans = function(s,s0) (s-s0)*log(s-s0)+exp(s-s0)
    y = scale(Ftrans(z,min(z)-1e-5))
    
    sel  =  which( abs( y ) > 10 )
    
    while (length( sel ) > 0) {
      y[ sel ] = NA
      y  =  scale( y )
      sel  =  which( abs( y ) > 10 )
    }
    
    idx = which(is.na(y))
    if (length(idx) != 0){
      y = y[-idx]
      GRS= as.matrix(GRS[-idx,])
    } else {
      y = y
      GRS = GRS
    }
    
    # Estimate interaction effect for GRS
    return( estimate_gxe( y, GRS, sim_num,simulate_phenotype = F) )
  })
  return(all)
}
stopCluster(cl)

allcomb5=allcomb
save(allcomb5,sim_xlogx_nor,file='./data/supp_xlogx_normal.Rdata')


######################### plot ################################

library(reshape2)
library(ggplot2)
library(gridExtra)

se = function(data){sd(data)/sqrt(length(data))}

rgamma = lapply(sim_exp,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][["gamma"]]))
rgamma = do.call(rbind,rgamma)
rgamma = cbind(data='GRS',as.data.frame(t(allcomb)),t(rgamma))
fgamma = lapply(sim_exp,function(x) sapply(x,function(n) n[["fake_grs"]][["coefficients"]][["gamma"]]))
fgamma = do.call(rbind,fgamma)
fgamma = cbind(data='FakeGRS',as.data.frame(t(allcomb)),t(fgamma))
gamma = rbind(rgamma,fgamma)

# long = melt(gamma,id.vars=c('data','c1','power','skwN','krtN'),variable.name = 'repeat',value.name = 'gamma')
# mean = aggregate(data=long,gamma~data+c1+power+krtN+skwN,'mean')
# sd = aggregate(data=long,gamma~data+c1+power+krtN+skwN,'se')
long = melt(gamma,id.vars=c('data','c1','skwN','krtN'),variable.name = 'repeat',value.name = 'gamma')
mean = aggregate(data=long,gamma~data+c1+krtN+skwN,'mean')
sd = aggregate(data=long,gamma~data+c1+krtN+skwN,'se')
try = cbind(mean,sd=sd$gamma)
colnames(try)[5] = 'mean'
try$data = factor(try$data)
try$c1 = factor(try$c1)
# try$power = factor(try$power)
try$skwN = factor(try$skwN)
## long[,1:4] = apply(long[,1:4],2,factor)

# try11 = try[which(try$c1==0 & try$power==-2),]
# try12 = try[which(try$c1==allc1[2] & try$power==-2),]
# try21 = try[which(try$c1==0 & try$power==-1),]
# try22 = try[which(try$c1==allc1[2] & try$power==-1),]

try11 = try[which(try$c1==0),]
try12 = try[which(try$c1==allc1[2]),]

plot= list()
plot[[1]] = ggplot(try11,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=0, f(t)=exp(t)')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[2]] = ggplot(try12,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=sqrt(0.05), f(t)=exp(t)')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[3]] = ggplot(try21,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=0, f(t)=t^-1')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[4]] = ggplot(try22,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=sqrt(0.05), f(t)=t^-1')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))

pdf('./results/supp_fig3.pdf',width=12,height=12)
grid.arrange(grobs=plot,nrow=1)
dev.off()





###################### supplimentary figure ##################
library(ggplot2)

load('./data/fig3.Rdata')   ### sim3
load('./data/supp_transform.Rdata') ### sim5
load('./data/supp_transform_exp.Rdata') ### sim_exp

rgamma = lapply(sim_exp,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][["gamma"]]))
rgamma = do.call(rbind,rgamma)
rgamma = cbind(data='GRS',as.data.frame(t(allcomb)),t(rgamma))
fgamma = lapply(sim_exp,function(x) sapply(x,function(n) n[["fake_grs"]][["coefficients"]][["gamma"]]))
fgamma = do.call(rbind,fgamma)
fgamma = cbind(data='FakeGRS',as.data.frame(t(allcomb)),t(fgamma))
gamma_exp = rbind(rgamma,fgamma)
gamma_exp = cbind(gamma_exp[,1:2],power='exp',gamma_exp[,-c(1,2)])

# save(sim3,gamma3,file='./data/fig3.Rdata')
# save(sim5,gamma5,file='./data/supp_transform.Rdata')
# save(sim_exp,gamma_exp,file='./data/supp_transform_exp.Rdata')

gamma = rbind(gamma3,gamma5,gamma_exp)
long = melt(gamma,id.vars=c('data','c1','power','skwN','krtN'),variable.name = 'repeat',value.name = 'gamma')
long[which(long$power=='0'),3] = 'log'
long$power = factor(long$power,levels = c('-2','-1','1','2','log','exp'))
longn = long[which(long$c1==0),]
longi = long[which(long$c1 != 0),]

supp = list()
supp[[1]] = ggplot(longn,aes(x=power,y=gamma,color=data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=0),linetype='dashed',color='grey')+
  labs(x='Transformation power (p) [f(t)=t^p]',y='estimates for gamma')
supp[[2]] = ggplot(longi,aes(x=power,y=gamma,color=data))+
  geom_boxplot()+
  geom_hline(aes(yintercept=sqrt(.05)),linetype='dashed',color='grey')+
  labs(x='Transformation power (p) [f(t)=t^p]',y='estimates for gamma')

pdf('./results/supp_all.pdf',width=12,height=4)
grid.arrange(grobs=supp,nrow=1)
dev.off()


############# significant p calculation #############
### when no transformation, diff between gamma & fake gamma
library(ggplot2)
library(reshape2)

## fig3 pars
allskwN = c(0,2,4)    # skewness of the noise
allkrtN = c(2,6,18) # kurtosis of the noise
allc1   = c(0,sqrt(0.05))        # gamma
allpow  = c(0,1,2) # transformation power
# allpow = 1

allcomb = expand.grid(allc1,allpow,allskwN,allkrtN)
allcomb = allcomb[which(allcomb[,4] > allcomb[,3]^2+1),]
colnames(allcomb) = c('c1','power','skwN','krtN')
allcomb = t(allcomb)

load('./data/fig3_0.05.Rdata')
load('./data/fig3_0.05_2.Rdata')

idx = which(allcomb[2,]==1)
pars = t(allcomb[,idx])

# notrans = lapply(sim32,function(x) x[idx])
# sim3 = c(notrans,sim32_2)
# p = lapply(sim3,function(x) sapply(x,function(n) n[["fake_grs"]][["p"]][["gamma"]]))
# p = do.call(rbind,p)
# sum = apply(p,2,function(x) sum(x<0.05)/length(x))
# mean(sum)
# res = cbind(pars,sum)
# res = res[order(res[,5],decreasing = T),]

### for gamma
aggregate(data=res,sum~c1,FUN = mean)

tidx = which(allcomb[2,]==0)
tpars = t(allcomb[,tidx])
trans = lapply(sim32,function(x) x[tidx])
tp = lapply(trans,function(x) sapply(x,function(n) n[["real_data"]][["p"]][["gamma"]]))
tp = do.call(rbind,tp)
tsum = apply(tp,2,function(x) sum(x<0.05)/length(x))
mean(tsum)
tres = cbind(tpars,tsum)


aggregate(data=tres,tsum~c1,FUN = mean) ### no sig diff

### t test

# fg = lapply(sim3,function(x) sapply(x,function(n) n[["fake_grs"]][["coefficients"]][["gamma"]]))
# fg = do.call(rbind,fg)
# t = apply(fg,2,function(x) t.test(x,mu=0)$p.value)
# mean = apply(fg,2,mean)

load('./data/fig3_0.05.Rdata')   ### sim3
load('./data/supp_transform.Rdata') ### sim5
load('./data/supp_transform_exp.Rdata') ### sim_exp

allgamma = list(gamma5,gamma3,gamma_exp)
alldiff = lapply(allgamma,function(gamma) gamma[which(gamma$data=='GRS'),6:55]-gamma[which(gamma$data!='GRS'),6:55])
alldiff = do.call(rbind,alldiff)
diff = apply(alldiff,1,function(x) sum(x>0)/length(x)) ##x>0: G>fG, larger percentage, more G>fG (what we need in real GxE)

allsim = list(pow2_1=sim5,pow0_2=sim3,exp=sim_exp)
allt = lapply(allsim,function(data) {
  t = lapply(data,function(x) sapply(x,function(n) n[["t_real_fgrs"]][["gamma"]]))
  t = do.call(rbind,t)
  return(t)
})
allt = do.call(cbind,allt)
allpar = rbind(gamma5[1:(nrow(gamma5)/2),2:5],gamma3[1:(nrow(gamma3)/2),2:5],gamma_exp[1:(nrow(gamma_exp)/2),2:5])
colnames(allpar)[2] = 'trans'
allpt = 2*pt(abs(allt),99,lower.tail = F)
power = apply(allpt,2,function(x) sum(x<0.05)/length(x))

allpard = cbind(allpar,diff)
wided = dcast(allpard,trans+skwN+krtN~c1)
colnames(wided)[4:5] = c('diff_0','diff_int')
allpart = cbind(allpar,power)
widet = dcast(allpart,trans+skwN+krtN~c1)
colnames(widet)[4:5] = c('power_0','power_int')
wide = cbind(wided,widet[,4:5])

## determine the symbol of power by difference of GxE
symbol = cbind(allpar,power,diff)
symbol$power[which(symbol$diff>=0.5)] = -symbol$power[which(symbol$diff>=0.5)]
wide2 = dcast(symbol,trans+skwN+krtN~c1,value.var='power')
colnames(wide2)[5] = 'sqrt(.05)'

# select normal skw and krt under GxE trans
comp_int = aggregate(data=allpart,power~c1+trans,mean)
comp_int$power = round(comp_int$power,2)
comp_int = dcast(comp_int,trans~c1)
colnames(comp_int)[3]= 'sqrt(.05)'
normal = widet[which(widet$skwN==0 & widet$krtN==2),]
similar = widet[which(widet$power_int<=0.5),]
diff = wide2[which(wide2$`sqrt(.05)`<= -0.5),]

## p of real grs gamma
allp = lapply(allsim,function(data) {
  p = lapply(data,function(x) sapply(x,function(n) n[["real_data"]][["p"]][["gamma"]]))
  p = do.call(rbind,p)
  return(p)
})
allp = do.call(cbind,allp)
perc_p = apply(allp,2,function(x) sum(x<0.05)/length(x))
allparp = cbind(symbol[1:5],perc_p)

sig = allparp[which(allparp$perc_p>0.5 & allparp$power< -0.5),]
sig = sig[order(sig$trans),]

##compare diff of beta & gamma of fY and real (c1=sqrt,skw=0,krt=2,pow=0)

fY_t = lapply(sim_pow1,function(x) {
  t_fy_beta = (x[["real_data"]][["coefficients"]][["beta"]] - x[["fake_phenotype"]][["coefficients"]][["beta"]]) / 
    sqrt(x[["real_data"]][["se"]][["beta"]]^2 + x[["fake_phenotype"]][["se"]][["beta"]])
  t_fy_gamma = (x[["real_data"]][["coefficients"]][["gamma"]] - x[["fake_phenotype"]][["coefficients"]][["gamma"]]) / 
    sqrt(x[["real_data"]][["se"]][["gamma"]]^2 + x[["fake_phenotype"]][["se"]][["gamma"]]^2)
  c(t_fy_beta,t_fy_gamma)
})
fY_t = do.call(cbind,fY_t)
fY_pt = 2*pt(abs(fY_t),99,lower.tail = F)
rowMeans(fY_pt)
perc_fY = apply(fY_pt,1,function(x) sum(x<0.05)/length(x))
### no significant in beta and gamma of fY --> so conclusion not satisfied

######### compare gamma of G and fG (with 0)
# load('./data/fig3_0.05.Rdata')   ### sim3
# load('./data/supp_transform.Rdata') ### sim5
# load('./data/supp_transform_exp.Rdata') ### sim_exp
# allsim = list(pow2_1=sim5,pow0_2=sim32,exp=sim_exp)

# allpar = rbind(gamma5[1:(nrow(gamma5)/2),2:5],allcomb,gamma_exp[1:(nrow(gamma_exp)/2),2:5]) # gamma3[1:(nrow(gamma3)/2),2:5]
# colnames(allpar)[2] = 'trans'

load('./data/fig3_0.05.Rdata')
load('./data/supp_pow1_normal.Rdata')
load('./data/supp_pow0_normal.Rdata')
load('./data/supp_pow2_normal.Rdata')
allsim = list(sim32,sim_1_nor,sim_0_nor,sim_2_nor)
allpar = rbind(allcomb,t(allcomb2),t(allcomb3),t(allcomb4))
colnames(allpar)[2] = 'trans'

allpg = lapply(allsim,function(data) {
  pg = lapply(data,function(x) sapply(x,function(n) n[["real_data"]][["p"]][["gamma"]]))
  pg = do.call(rbind,pg)
  return(pg)
})
allpg = do.call(cbind,allpg)
allpfg = lapply(allsim,function(data) {
  pfg = lapply(data,function(x) sapply(x,function(n) n[["fake_grs"]][["p"]][["gamma"]]))
  pfg = do.call(rbind,pfg)
  return(pfg)
})
allpfg = do.call(cbind,allpfg)

perc_sig = colSums(allpg< 0.05 & allpfg>= 0.05)/nrow(allpg)
perc_pg = apply(allpg,2,function(x) sum(x<0.05)/length(x))
perc_pfg = apply(allpfg,2,function(x) sum(x<0.05)/length(x))
res = cbind(allpar,perc_pg,perc_pfg,perc_sig)
res[which(res$trans !=1 & res$perc_sig>0.5),] ##situation 2

allt = lapply(allsim,function(data) {
  t = lapply(data,function(x) sapply(x,function(n) n[["t_real_fgrs"]][["gamma"]]))
  t = do.call(rbind,t)
  return(t)
})
allt = do.call(cbind,allt)
allpt = 2*pt(abs(allt),99,lower.tail = F)
perc_pt = apply(allpt,2,function(x) sum(x>=0.05)/length(x))
perc_sig2 = colSums(allpg<0.05 & allpt>=0.05)/nrow(allpg)
res2 = cbind(allpar,perc_pg,perc_pt,perc_sig2)
res2[which(res2$trans==1 & res2$perc_sig2>0.5),]


## fig3 pars
allskwN = c(0,2,4)    # skewness of the noise
allkrtN = c(2,6,18) # kurtosis of the noise
allc1   = c(0,sqrt(0.05))        # gamma
allpow  = c(0,1,2) # transformation power
# allpow = 1

allcomb = expand.grid(allc1,allpow,allskwN,allkrtN)
allcomb = allcomb[which(allcomb[,4] > allcomb[,3]^2+1),]
colnames(allcomb) = c('c1','power','skwN','krtN')

pg = lapply(sim32,function(x) sapply(x,function(n) n[["real_data"]][["p"]][["gamma"]]))
pg = do.call(rbind,pg)
pfg = lapply(sim32,function(x) sapply(x,function(n) n[["fake_grs"]][["p"]][["gamma"]]))
pfg = do.call(rbind,pfg)
##method 1
perc_sig = colSums(pg< (0.05) & pfg>= (0.05))/nrow(pg)
res = cbind(allcomb,perc_sig)
res[which(res$perc_sig>0.5),]
##method 2
res = cbind(allcomb,pg=colMeans(pg),pfg=colMeans(pfg))
res[which(res$pg<0.05 & res$pfg>=0.05),]
##method 3
perc_pg = apply(pg,2,function(x) sum(x<0.05)/length(x))
perc_pfg = apply(pfg,2,function(x) sum(x<0.05)/length(x))
res = cbind(allcomb,perc_pg,perc_pfg)
sig = res[which(res$perc_pg>0.5 & res$perc_pfg<0.5),]

res = cbind(allcomb,perc_pg,perc_pfg,perc_sig)
res[which(res$power !=1 & res$perc_sig>0.5),]

# tidx = which(allcomb[2,]!=100)
# tpars = t(allcomb[,tidx])
# 
# data = lapply(sim32,function(x) x[tidx])
# data = sim3
# 
# t = lapply(data,function(x) sapply(x,function(n) n[["t_real_fgrs"]][["gamma"]]))
# t = do.call(rbind,t)
# pt = 2*pt(abs(t),99,lower.tail = F)
# diff = apply(pt,2,function(x) sum(x<0.05)/length(x))
# mean(diff)
# part = cbind(tpars,diff)
# part = as.data.frame(part)
# # part = part[order(part$power,part$c1),]
# try = dcast(part,power+skwN+krtN~c1)
# colnames(try)[5] = 'sqrt(.05)'
# 
# 
# compare = aggregate(data=part,diff~c1+power,'mean')
# compare2 = aggregate(data=part,diff~power,'mean')
# 
# 
# part$c1 = factor(part$c1)
# part$skwN = factor(part$skwN)
# part$krtN = factor(part$krtN)
# ggplot(part,aes(x=c1,y=diff,color=krtN,shape=skwN))+geom_point()

# # sum(t<0.05)/length(t)
# # fg2 = as.numeric(fg)
# # t.test(fg2,mu=0)$p.value
# ggplot(part,aes(x=krtN,y=diff,color=skwN,shape=c1))+geom_point()
