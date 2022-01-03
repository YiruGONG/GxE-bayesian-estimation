computer = 0    ### 1 for own computer ; 0 for thinkpad
dirct = ifelse(computer,'D:/Jiang Lab/','F:/')
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


############ Fig 3 - Outcome modelled on transformed scale ############
## n = 10000, a1^2=0.1,a2=0, beta^2=0.3, E~N(0,1),f(t)=t
##  gamma^2=0 or gamma^2=0.05
## transformation f(t)=t^p, p=0,1,2
## generate fake GRS

m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.05)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E

allskwN = c(0,2,4)    # skewness of the noise
allkrtN = c(2,6,18) # kurtosis of the noise
allc1   = c(0,sqrt(0.05))        # gamma
allpow  = c(0,1,2) # transformation power
# allpow = 1

allcomb = expand.grid(allc1,allpow,allskwN,allkrtN)
allcomb = allcomb[which(allcomb[,4] > allcomb[,3]^2+1),]
colnames(allcomb) = c('c1','power','skwN','krtN')
allcomb = t(allcomb)

# Simulation
###open multiple cores
cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)

sim32_2 <- foreach(rep=1:100,.export=c('rpearson','estimate_gxe','mclapply')) %dopar% {
  all <- lapply(as.data.frame(allcomb),function(x){
    c1 = x[1]
    pow= x[2]
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
    
    Ftrans  =  function( s, p1, p2 ) ( (s-p1)^p2 - 1 ) / p2 
    
    if (pow != 0) {
      y  =  scale( Ftrans( z, min(z)-1e-5, pow ) )
    } else {
      y  =  scale( log( z - min(z)+1e-5 ) )
    }
    
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
    return( estimate_gxe( y, GRS, sim_num) )
  })
  return(all)
}
stopCluster(cl)

save(sim32_2,file='./data/fig3_0.05_2.Rdata')

library(reshape2)
library(ggplot2)
library(gridExtra)

se = function(data){sd(data)/sqrt(length(data))}

rgamma = lapply(sim32,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][["gamma"]]))
rgamma = do.call(rbind,rgamma)
rgamma = cbind(data='GRS',as.data.frame(t(allcomb)),t(rgamma))
fgamma = lapply(sim3,function(x) sapply(x,function(n) n[["fake_grs"]][["coefficients"]][["gamma"]]))
fgamma = do.call(rbind,fgamma)
fgamma = cbind(data='FakeGRS',as.data.frame(t(allcomb)),t(fgamma))
gamma = rbind(rgamma,fgamma)

long = melt(gamma,id.vars=c('data','c1','power','skwN','krtN'),variable.name = 'repeat',value.name = 'gamma')
mean = aggregate(data=long,gamma~data+c1+power+krtN+skwN,'mean')
sd = aggregate(data=long,gamma~data+c1+power+krtN+skwN,'se')
try = cbind(mean,sd=sd$gamma)
colnames(try)[6] = 'mean'
try$data = factor(try$data)
try$c1 = factor(try$c1)
try$power = factor(try$power)
try$skwN = factor(try$skwN)
# long[,1:4] = apply(long[,1:4],2,factor)

try11 = try[which(try$c1==0 & try$power==0),]
try12 = try[which(try$c1==allc1[2] & try$power==0),]
try21 = try[which(try$c1==0 & try$power==1),]
try22 = try[which(try$c1==allc1[2] & try$power==1),]
try31 = try[which(try$c1==0 & try$power==2),]
try32 = try[which(try$c1==allc1[2] & try$power==2),]

plot= list()
plot[[1]] = ggplot(try11,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=0, f(t)=log(t)')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[2]] = ggplot(try12,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=sqrt(0.05), f(t)=log(t)')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[3]] = ggplot(try21,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=0, f(t)=t')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[4]] = ggplot(try22,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=sqrt(0.05), f(t)=t')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[5]] = ggplot(try31,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=0, f(t)=t^2')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))
plot[[6]] = ggplot(try32,aes(krtN,mean,color=skwN,shape=data))+geom_point(size=3)+ #,position=position_dodge(0.2))+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,size=0.75)+ #position=position_dodge(0.2))+
  # scale_y_continuous(breaks=seq(-0.3,0,0.02))+
  scale_x_continuous(breaks=allkrtN)+
  labs(x='Kurtosis',y='estimated gamma',title='gamma=sqrt(0.05), f(t)=t^2')+
  theme_bw()+
  scale_color_manual(values=c('blue','green','red'))

pdf('./results/fig32.pdf',width=12,height=12)
grid.arrange(grobs=plot,nrow=3)
dev.off()
