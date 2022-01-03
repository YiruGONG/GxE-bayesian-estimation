computer = 1    ### 1 for own computer ; 2 for thinkpad
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

################ Fig 4 - power analysis ######################
## gamma^2 =seq(0.002,0.02,0.002), sample size (m) = seq(1e4,10e4,1e4)

m    = 100        # number of genetic markers
a1   = sqrt(.1)   # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)   # beta
skwE = 0          # skewness of E
krtE = 3          # kurtosis of E
skwN = 0          # skewness of the noise
krtN = 3          # kurtosis of the noise
pow  = 1          # transformation power

alln    = seq(1e4,1e5,2e4) # sample size
allc1   = sqrt(seq(0.002,0.02,0.004)) # gamma

allcomb = expand.grid(alln,allc1)
colnames(allcomb) = c('sample','gamma')
allcomb = t(allcomb)

# Simulation
###open multiple cores
cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)

sim4 <- foreach(rep=1:50,.export=c('rpearson','estimate_gxe','mclapply')) %dopar% {
  # start=Sys.time()
  all <- lapply(as.data.frame(allcomb),function(x){
    n = as.numeric(x[1])
    c1 = as.numeric(x[2])
    
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
      GRS = GRS}
    
    # Estimate interaction effect for GRS
    return( estimate_gxe( y, GRS, sim_num) )
  })
  # Sys.time()-start
  save(all,file=paste0('./data/fig4/rep_',rep,'.Rdata'))
  return(c(all,n=n,c1=c1))
}
stopCluster(cl)

save(sim4,file='./data/fig4.Rdata')

# se = function(data){sd(data)/sqrt(length(data))}
library(ggplot2)
library(RColorBrewer)

col = scales::seq_gradient_pal('darkgreen','yellow')(seq(0,1,length.out = 5))

load('./data/fig4.Rdata')
pgamma = lapply(sim4,function(x) sapply(x,function(n) n[["real_data"]][["p"]][["gamma"]]))
pgamma = do.call(rbind,pgamma)
power = apply(pgamma,2,function(x) sum(x<1e-3)/length(x))
se = lapply(sim4,function(x) sapply(x,function(n) n[["real_data"]][["se"]][["gamma"]]))
se = do.call(rbind,se)
se = apply(se,2,mean)
power = cbind(t(allcomb),as.data.frame(power),se)
power$sample = factor(power$sample)
power$gamma2 = power$gamma^2

pdf('./results/fig4.pdf',width=5,height = 4)
ggplot(power,aes(gamma2,power,color=sample))+geom_line()+geom_point(size=2)+
  geom_errorbar(aes(ymin=power-se,ymax=power+se,color=sample),width=.0003)+
  scale_color_manual(values=col)+
  labs(x='gamma^2',y='Power (at 1e-3)')+
  scale_x_continuous(breaks = allc1^2)
dev.off()

# par= 'beta'
# se = aa[["real_data"]][["se"]][par]
# coef_idv = aa[["real_data"]][["individual_estimates"]][par,]
# p_idv =  2 * pnorm( -abs(  coef_idv / se ) )
# power = sum(p_idv < 1e-3) / length(p_idv)

