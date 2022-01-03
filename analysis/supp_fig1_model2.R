setwd('D:/Jiang Lab/GxE')

# library('devtools')
# install_github('zkutalik/GRSxE_software',subdir='Rcode')
# library( GxE )
library( PearsonDS )
library(doParallel)
library(foreach)

source('D:/Jiang Lab/GxE/GRSxE_software/Rcode/R/estimate_gxe.R')

############## Fig 1 - Effect of G¨CE correlation on the parameter estimation ##############
##correlation range from 1~0.3,seq=0.05
##500 simulated data sets
## n = 10000, alpha'=0.1, beta'^2=0.3, gamma'^2=0.05, E~N(0,sig^2)

m    = 100       # number of genetic markers
n    = 1e4       # sample size
a1   = sqrt(.1)  # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)  # beta
# c1   = sqrt(.3)  # gamma
skwE = 0         # skewness of E
krtE = 3         # kurtosis of E
skwN = 0         # skewness of the noise
krtN = 3         # kurtosis of the noise
pow  = 1         # transformation power

# b   = seq(0,0.4,0.1)
c   = sqrt( seq(0,0.4,0.1) )        # correlation between E and GRS

# d=0.3

if (krtE <= skwE^2 + 1 | krtN <= skwN^2 + 1) {
  stop( 'Skew and kurtosis values not compatible (kurtosis > skew^2 + 1 not satisfied)' )
}

sim_num = 1e2 # number of bootstrap / fake GRS

# Simulation



cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)
# clusterExport(cl,c(),envir = environment())

# sim1_0.5 = parLapply(cl,1:100,function(x){
sim1_m2_c = foreach(rep=1:100,.export=c('rpearson','estimate_gxe','mclapply')) %dopar% {
  # start=Sys.time()
  allc = lapply(t(as.matrix(c)),function(c1){
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
    return( estimate_gxe( y, GRS, sim_num ) )
  })
  # Sys.time()-start
  return(allc)
}
stopCluster(cl)
Sys.time()
save(sim1_m2_c,file='./data/fig1_m2_c.Rdata')



beta = lapply(sim1_0.5,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][["beta"]]))
beta = do.call(rbind,beta)
colnames(beta) = d
gamma = lapply(sim1_0.5,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][["gamma"]]))
gamma = do.call(rbind,gamma)
colnames(gamma) = d
delta = lapply(sim1_0.5,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][["alpha2"]]/c1))
delta = do.call(rbind,delta)
colnames(delta) = d

par1 = lapply(sim1_0.5,function(x) lapply(x,function(n) {
  b=n[["real_data"]][["coefficients"]][["alpha2"]]
  c=n[["real_data"]][["coefficients"]][["beta"]]
  d=n[["real_data"]][["coefficients"]][["gamma"]]
  gamma0= sqrt(b^2+d^2)
  delta = b/gamma0
  beta0 = c/sqrt(1-delta^2)
  return(c(beta0=beta0,gamma0=gamma0,delta=delta))}))
beta0 = lapply(par1,function(x) sapply(x,function(n) n[["beta0"]]))
beta0 = do.call(rbind,beta0)
gamma0 = lapply(par1,function(x) sapply(x,function(n) n[["gamma0"]]))
gamma0 = do.call(rbind,gamma0)
delta = lapply(par1,function(x) sapply(x,function(n) n[["delta"]]))
delta = do.call(rbind,delta)
colnames(beta0) <- colnames(gamma0) <- colnames(delta) <- d

pdf('./results/fig1_0.5_2.pdf',width = 5,height=4)
boxplot(beta0,xlab='delta',ylab='Estimation for beta',las=1)
abline(h=b1,lty=2)
boxplot(gamma0,xlab='delta',ylab='Estimation for gamma',las=1)
abline(h=c1,lty=2)
boxplot(delta,xlab='delta',ylab='Estimation for delta',yaxt='n')
axis(2,at=seq(-0.1,0.4,0.05),las=1)
arrows(c(1:7)-0.6,d,c(1:7)+0.6,d,angle=0,lty=2)
# abline(h=d,lty=2,col='grey')
dev.off()

# ##test for bootsrap
# try = lapply(sim1_0.5,function(x) sapply(x,function(n) n[["real_data"]][["individual_estimates"]]['beta',]))
# try = do.call(rbind,try[1:5])
# colnames(try) = d
# boxplot(try)
# abline(h=b1,lty=2)
