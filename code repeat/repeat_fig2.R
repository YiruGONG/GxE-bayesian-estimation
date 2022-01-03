computer = 1    ### 1 for own computer ; 2 for thinkpad
dirct = ifelse(computer,'D:/Jiang Lab/','F:/')
setwd(paste0(dirct,'GxE'))

# library( GxE )
library( PearsonDS )
library(doParallel)
library(foreach)

##################### Effect of non-normality of Error ###################

### the effect of skewness and kurtosis of the environmental variable (E) on the estimation bias
## n = 10000, a1^2=0.1,a2=0, beta'^2=0.3, gamma'^2=0, E~N(0,1),f(t)=t
## diff E(e^3) and E(e^4)

m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
c1   = 0         # gamma
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E
allskwN = c(0,2,4) # skewness of the noise
allkrtN = c(2,6,18) # kurtosis of the noise
pow  = 1         # transformation power

allcomb = expand.grid(allskwN,allkrtN)
allcomb = allcomb[which(allcomb[,2] > allcomb[,1]^2+1),]
colnames(allcomb) = c('skwN','krtN')
allcomb = t(allcomb)


# Simulation
###open multiple cores
cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)

sim21 <- foreach(comb=1:ncol(allcomb),.export=c('rpearson','estimate_gxe','mclapply')) %dopar% {
  # sim21 = parLapply(cl,allcomb,function(x){
  skwN = allcomb[1,comb]
  krtN = allcomb[2,comb]
  
  all <- lapply(1:100,function(x){
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
    if (length(idx) != 0) {
      y= y[-idx]
      GRS = as.matrix(GRS[-idx,])
    } else {
      y=y
      GRS=GRS}
    

    
    # Estimate interaction effect for GRS
    return( estimate_gxe( y, GRS, sim_num ) )
  })
  return(c(all,skwN,krtN))
}
stopCluster(cl)




### effect of skewness and kurtosis of the residual noise (ϵ) on the estimation bias 
## n = 10000, a1^2=0.1,a2=0, beta'^2=0.3, gamma'^2=0, e~N(0,1),f(t)=t
## diff E(e^3) and E(e^4)

m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
c1   = 0         # gamma
allskwE = c(0,2,4)    # skewness of E
allkrtE = c(2,6,18) # kurtosis of E
skwN = 0         # skewness of the noise
krtN = 3         # kurtosis of the noise
pow  = 1         # transformation power

allcomb = expand.grid(allskwE,allkrtE)
allcomb = allcomb[which(allcomb[,2] > allcomb[,1]^2+1),]
colnames(allcomb) = c('skwE','krtE')
allcomb = t(allcomb)

# Simulation
###open multiple cores
cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)

sim22 <- foreach(comb=1:ncol(allcomb),.export=c('rpearson','estimate_gxe','mclapply')) %dopar% {
  # sim21 = parLapply(cl,allcomb,function(x){
  skwE = allcomb[1,comb]
  krtE = allcomb[2,comb]
  
  all <- lapply(1:50,function(x){
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
    if (length(idx) != 0) {
      y= y[-idx]
      GRS = as.matrix(GRS[-idx,])
    } else {
      y=y
      GRS=GRS}
    
    # Estimate interaction effect for GRS
    return( estimate_gxe( y, GRS, sim_num ) )
  })
  return(all)
}
stopCluster(cl)
save(sim21,sim22,file='./data/fig2.Rdata')

# sim=list()
# for (comb in 1:length(sim22)){
#   rep = c(sim22[[comb]],sim23[[comb]])
#   sim[[comb]] = c(rep,par=allcomb[,comb]) 
# }
# sim22=sim

library(reshape2)
library(ggplot2)
library(gridExtra)

BE = lapply(sim22,function(x) c(x$par.skwN,x$par.krtN, sapply(x[1:100], function(n) n[["real_data"]][["coefficients"]][["gamma"]])))
BE = do.call(rbind,BE)
BE = as.data.frame(BE)
BE[,1] = factor(BE[,1])
BE[,2] = factor(BE[,2],allkrtE)
colnames(BE) = c('skwE','krtE',1:100)
Bn = lapply(sim21,function(x) c(x[[101]],x[[102]], sapply(x[1:100], function(n) n[["real_data"]][["coefficients"]][["gamma"]])))
Bn = do.call(rbind,Bn)
Bn = as.data.frame(Bn)
Bn[,1] = factor(Bn[,1])
Bn[,2] = factor(Bn[,2],allkrtN)
colnames(Bn) = c('skwN','krtN',1:100)

BEs = melt(BE[,-2],id.vars=c('skwE'),varnames ='rep',value.name = 'gamma')
BEk = melt(BE[,-1],id.vars=c('krtE'),varnames ='rep',value.name = 'gamma')
Bns = melt(Bn[,-2],id.vars=c('skwN'),varnames ='rep',value.name = 'gamma')
Bnk = melt(Bn[,-1],id.vars=c('krtN'),varnames ='rep',value.name = 'gamma')

plot=list()
plot[[1]]= ggplot(BEs,aes(skwE,gamma))+geom_boxplot()+
  geom_hline(yintercept = 0,linetype='dashed',color='red')+
  labs(x='Skewness(E)',y='Bias')
plot[[2]]= ggplot(BEk,aes(krtE,gamma))+geom_boxplot()+
  geom_hline(yintercept = 0,linetype='dashed',color='red')+
  labs(x='Kurtosis(E)',y='Bias')
plot[[3]] = ggplot(Bns,aes(skwN,gamma))+geom_boxplot()+
  geom_hline(yintercept = 0,linetype='dashed',color='red')+
  labs(x='Skewness(e)',y='Bias')
plot[[4]]= ggplot(Bnk,aes(krtN,gamma))+geom_boxplot()+
  geom_hline(yintercept = 0,linetype='dashed',color='red')+
  labs(x='Kurtosis(e)',y='Bias')

pdf('./results/fig2.pdf',width = 8,height=6)
grid.arrange(grobs=plot,nrow=2)
dev.off()

