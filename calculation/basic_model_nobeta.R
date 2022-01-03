setwd('D:/Jiang Lab/GxE')

library( PearsonDS )
# library(RNOmni)
library(doParallel)
library(foreach)
# library(gsl)

source('./calculation/basic_fit.R')
source('./calculation/basic_fit_G2.R')

# source( './GRSxE_software/Rcode/R/IA_fit.R' )
# source( './GRSxE_software/Rcode/R/IA_fit_G2.R' )


.single_gxe  =  function( y,
                          grs,
                          params ){
  params  =  optim( params,
                    basic_fit,
                    y = y, grs = grs )$par
  optim( c( params[ 1 ], 0, params[ 2:4 ] ),
         basic_fit_G2,
         y = y, grs = grs )$par
}

estimate_basic = function(y,grs){
  model = lm( y ~ grs + I( grs^2) )
  init = model$coefficients[ 2 ]
  var = sum((model$residuals)^2)/(n-3) 
  init_var = sqrt(var/2)
  coeff = .single_gxe(y,grs,c(init,0,init_var,init_var))
  coef_names = c('alpha1', 'alpha2', 'gamma','sigmaE','sigmaN' )
  return(setNames( coeff, coef_names ))
}



################ simulation ###################

m    = 100       # number of genetic markers
n    = 1e4       # sample size
# a0   = 0.1
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
# sigE = sqrt(.5)
c1   = sqrt(.3)  # gamma
sigN = sqrt(.5)
sigE = sqrt(.5)
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E
skwN = 0         # skewness of the noise
krtN = 3         # kurtosis of the noise
# pow  = 1

b   = sqrt(seq(0,0.3,0.05))  # beta
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

basic_nobeta = foreach(b1=b,.export=c('rpearson')) %dopar% {
  # start=Sys.time()
  # alld = lapply(t(as.matrix(d)),function(d1){
  reps = sapply(1:500,function(rep){
    
    maf  =  matrix( runif( m ), nrow = 1 )
    G    =  apply( maf, 2, function(x) rbinom( n, 2, x ) )
    eff  =  matrix( runif( m ), ncol = 1 )
    eff  =  eff / sqrt( sum( eff^2 ) )
    GRS  =  scale( G %*% eff )
    
    E    =  matrix( rpearson( n, moments = c( 0, 1, skwE, krtE ) ) )
    noi  =  matrix( rpearson( n, moments = c( 0, 1, skwN, krtN ) ) )
    # sig  =  sqrt( 1 - a1^2 - 2*a2^2 - b1^2 - c1^2 )   ###transformation of sig
    
    # z    =  scale( a1 * GRS + a2*(GRS^2-1) + b1 * E + c1 * GRS * E + sig * noi )
    z    =  a1 * GRS + a2*(GRS^2-1) + b1 *sigE* E + c1 * GRS * sigE*E + sigN * noi 
    
    y = z - mean(z)
    
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
    return( estimate_basic( y, GRS ) )
  })
  # Sys.time()-start
  return(reps)
}
stopCluster(cl)
Sys.time()
save(basic_nobeta,file='./data/new/basic_nobeta.Rdata')

# basic2 = lapply(basic,function(data) {
#   data = as.data.frame(t(data))
#   gamma0= sqrt(data$alpha2^2+data$gamma^2)
#   delta = data$alpha2/gamma0
#   beta0 = data$beta/sqrt(1-delta^2)
#   alpha0= data$alpha1-beta0*delta
#   return(cbind(alpha0,beta0,gamma0,delta))
# })
# 
# alpha0 = sapply(basic2,function(x) x[,'alpha0'])
# beta0 = sapply(basic2,function(x) x[,'beta0'])
# gamma0 = sapply(basic2,function(x) x[,'gamma0'])
# delta = sapply(basic2,function(x) x[,'delta'])
# colnames(alpha0) <- colnames(beta0) <- colnames(gamma0) <- colnames(delta) <- d

# alpha0 = sapply(basic_noscale,function(x) x['alpha0',])
alpha1 = sapply(basic_nobeta,function(x) x['alpha1',])
alpha2 = sapply(basic_nobeta,function(x) x['alpha2',])
# beta = sapply(basic_noscale,function(x) x['beta',])
gamma = sapply(basic_nobeta,function(x) x['gamma',])
sigmaE = sapply(basic_nobeta,function(x) x['sigmaE',])
sigmaN = sapply(basic_nobeta,function(x) x['sigmaN',])
colnames(alpha1) <- colnames(alpha2) <- colnames(sigmaE) <- colnames(sigmaN) <- colnames(gamma) <- b^2

pdf('./results/new/basic_nobeta.pdf',width = 10,height=12)
par(mfrow=c(3,2))
# boxplot(alpha0,xlab='different beta^2',ylab='Estimation for a0',las=1)
# abline(h=a0,lty=2)
boxplot(alpha1,xlab='different beta^2',ylab='Estimation for alpha1',main='alpha1',las=1)
abline(h=a1,lty=2)
boxplot(alpha2,xlab='different beta^2',ylab='Estimation for alpha2',main='alpha2',las=1)
abline(h=a2,lty=2)
# boxplot(beta,xlab='different beta^2',ylab='Estimation for beta',las=1)
# arrows(c(1:7)-0.6,b,c(1:7)+0.6,b,angle=0,lty=2)
boxplot(gamma,xlab='different beta^2',ylab='Estimation for gamma',main='gamma',las=1,ylim=c(0,4))
abline(h=c1,lty=2)
boxplot(sigmaE,xlab='different beta^2',ylab='Estimation for sigma E',main='sigma E',las=1,ylim=c(0,1))
abline(h=sigE,lty=2)
boxplot(sigmaN,xlab='different beta^2',ylab='Estimation for sigma N',main='sigma N',las=1,outline=F)
abline(h=sigN,lty=2)
boxplot(gamma,xlab='different beta^2',ylab='Estimation for gamma',main='Supplementary gamma',las=1)
abline(h=c1,lty=2)
dev.off()

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
c1   = sqrt(.3)  # gamma
# skwE = 0         # skewness of E
# krtE = 3         # kurtosis of E
# allskwN = c(0,2,4) # skewness of the noise
# allkrtN = c(2,6,8,18) # kurtosis of the noise
pow  = 1         # transformation power

skwN = 0         # skewness of E
krtN = 3         # kurtosis of E
allskwE = c(0,2,4) # skewness of the noise
allkrtE = c(2,6,8,18) # kurtosis of the noise

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
basic_skwE <- foreach(comb=1:ncol(allcomb),.export=c('rpearson')) %dopar% {
  # sim21 = parLapply(cl,allcomb,function(x){
  skwE = allcomb[1,comb]
  krtE = allcomb[2,comb]
  
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
    return( estimate_basic( y, GRS ) )
  })
  # Sys.time()-start
  return(reps)
}
stopCluster(cl)
Sys.time()
save(basic_skwE,file='./data/new/basic_skwE.Rdata')

# a0 = sapply(basic_skwE,function(x) x['a0',])
alpha1 = sapply(basic_skwE,function(x) x['alpha1',])
# alpha1 = cbind(t(allcomb),alpha1)
alpha2 = sapply(basic_skwE,function(x) x['alpha2',])
beta = sapply(basic_skwE,function(x) x['beta',])
gamma = sapply(basic_skwE,function(x) x['gamma',])
# sigE = sapply(basic_skwE,function(x) x['sigmaE',])

pdf('./results/new/basic_skwE.pdf',width = 10,height=8)
par(mfrow=c(2,2))
boxplot(alpha1,xlab='skw and krt',ylab='Estimation for alpha1',las=1)
abline(h=a1,lty=2)
boxplot(alpha2,xlab='skw and krt',ylab='Estimation for alpha2',las=1)
abline(h=a2,lty=2)
boxplot(beta,xlab='skw and krt',ylab='Estimation for beta',las=1)
abline(h=b1,lty=2)
boxplot(gamma,xlab='skw and krt',ylab='Estimation for gamma',las=1)
abline(h=c1,lty=2)
dev.off()