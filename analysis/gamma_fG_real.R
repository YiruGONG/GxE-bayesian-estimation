setwd('D:/Jiang Lab/GxE')

library( PearsonDS )
library(RNOmni)
library(doParallel)
library(foreach)
library(gsl)

source('./GRSxE_software/Rcode/R/estimate_gxe.R')



################ supp - compare real data gamma~e3+e4 #############

m    = 100        # number of genetic markers
n    = 1e4
a1   = sqrt(.1)   # linear effect of GRS on y
a2   = 0         
b1   = sqrt(.3)   # beta
skwE = 0          # skewness of E
krtE = 3          # kurtosis of E
skwN = 0          # skewness of the noise
krtN = 3          # kurtosis of the noise
pow  = 1          # transformation power

allc1   = sqrt(seq(0.002,0.02,0.002)) # gamma

# Simulation
###open multiple cores
cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(cl)

supp_r <- foreach(rep=1:50,.export=c('rpearson','estimate_gxe','mclapply')) %dopar% {
  # start=Sys.time()
  all <- lapply(allc1,function(x){
    c1 = x
    
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
    
    ### 4 es ###
    e3 = mean(GRS*y)
    e4 = mean(GRS^2*y)
    
    # Estimate interaction effect for GRS
    return( c(estimate_gxe( y, GRS, sim_num),list(c(e3,e4)) ) )
  })
  # Sys.time()-start
  # save(all,file=paste0('./data/fig4/rep_',rep,'.Rdata'))
  return(all)
}
stopCluster(cl)

save(supp_r,file='./data/supp_gamma_fG_real.Rdata')

###### compare ###########
e3s = lapply(supp_r,function(x) sapply(x,function(n) n[[5]][1]))
e4s = lapply(supp_r,function(x) sapply(x,function(n) n[[5]][2]))
gammas = lapply(supp_r,function(x) sapply(x,function(n) n[["real_data"]][["coefficients"]][["gamma"]]))
e3s = do.call(rbind,e3s)
e4s = do.call(rbind,e4s)
gammas = do.call(rbind,gammas)

# e3 = apply(e3s,2,mean)
# e4 = apply(e4s,2,mean)
gamma = apply(gammas,2,mean)
e3 = as.numeric(e3s)
e4 = as.numeric(e4s)
gamma = as.numeric(gammas)


var = outer(cbind(e3,e4),c(1,2,3,4),function(x,y) x^y)
size = dim(var)
dim(var) = c(size[1],size[2]*size[3])
colnames(var) = c('e3.1','e4.1','e3.2','e4.2','e3.3','e4.3','e3.4','e4.4')
int = expand.grid(c('e3.1','e3.2','e3.3','e3.4'),c('e4.1','e4.2','e4.3','e4.4'))
##??? t(int)
var_int = apply(t(int),2,function(x) var[,x[1]]*var[,x[2]]) ##???
colnames(var_int) = paste(int[,1],int[,2],sep='_')
var = cbind(gamma,var,var_int)
var = as.data.frame(var)

res1 = lm(gamma~.,data=var)
summary(res1)
res.bic = step(res1,k=log(nrow(var)),trace = F)
summary(res.bic)


