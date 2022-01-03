# computer = 1    ### 1 for own computer ; 0 for thinkpad
# dirct = ifelse(computer,'D:/Jiang Lab/','F:/')
# setwd(paste0(dirct,'GxE'))
setwd('D:/Jiang Lab/GxE')

library( PearsonDS )
library(RNOmni)
library(doParallel)
library(foreach)
library(gsl)
library( parallel )

source( 'D:/Jiang Lab/GxE/GRSxE_software/Rcode/R/IA_fit.R' )
source( 'D:/Jiang Lab/GxE/GRSxE_software/Rcode/R/IA_fit_G2.R' )


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

allc1   = c(0,sqrt(0.05))        # gamma
range=seq(-0.25,0.25,0.05)

results = list()
for (nn in 1){
  c1 = allc1[nn]
  
  ###### simulation #######
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
  
  ####### Estimate interaction effect for GRS ########
  # return( estimate_gxe( y, GRS, sim_num ) )
  
  .single_gxe  =  function( y,
                            grs,
                            params ){
    params  =  optim( params,
                      IA_fit,
                      y = y, grs = grs )$par
    optim( c( params[ 1 ], 0, params[ 2:3 ] ),
           IA_fit_G2,
           y = y, grs = grs )$par
  }
  
########### generate grs and y ##################
  
  # gamma_fG  =  function( y,
  #                            grs,
  #                            sim_num = 100,
  #                            simulate_phenotype = FALSE,
  #                            # The rest is only used if simulate_phenotype is TRUE
  #                            skewness_range = seq( -3, 3, by = 0.2 ),
  #                            k_range = c( 2:4 ), # kurtosis = skewness^2 + k
  #                            max_sd = 7 ) {
    if (is.null( dim( y ) )) {
      y = as.matrix( y )
    }
    
  grs = GRS
    ao_het  =  lm( y ~ grs + I( grs^2 ) )$coefficients[ 2 ]  ###starting point
    individual_coefficients  =  mclapply( 1:sim_num,
                                          function(x) {
                                            ix  =  sample( length(y),
                                                           length(y),
                                                           replace = TRUE )
                                            .single_gxe( y[   ix, , drop = FALSE ],
                                                         grs[ ix, , drop = FALSE ],
                                                         c( ao_het, 0.1, 0 ) )
                                          } )
    coef_names = c( 'alpha1', 'alpha2', 'beta', 'gamma' )
    individual_coefficients      =  do.call( cbind, individual_coefficients  )
    rownames( individual_coefficients      )  =  coef_names
    coefficients         =  apply( individual_coefficients,      1, mean )
    se_coefficients      =  apply( individual_coefficients,      1, sd   )
    p_coefficients       =  2 * pnorm( -abs( coefficients      / se_coefficients      ) )
    
    ######### generate fake GRSs ###############
    # sim_fGRS = function(y,grs,sim_num,range=seq(-0.25,0.25,0.05)){
      ge1 = mean(grs)
      ge2 = mean(grs^2)
      ge3 = mean(grs*y)
      ge4 = mean(grs^2*y)
      
      e1 = 0
      e2 = 1
      e3s = ge3 + range
      e4s = ge4 + range
      allcomb = expand.grid(e3s,e4s)
      
      ###open multiple cores
      cl.cores = detectCores(logical = F)
      cl <- makeCluster(cl.cores)
      registerDoParallel(cl)
      
      fGRS <- foreach(comb=1:nrow(allcomb),.export = '.single_gxe') %dopar% {
        e3 = allcomb[comb,1]
        e4 = allcomb[comb,2]
        
        mu1 = 0
        mu2 = 1
        mu3 = mean(y^3)
        mu4 = mean(y^4)
        mu5 = mean(y^5)
        
        A = mu3^2 + mu5 -2*mu3*mu4
        B = e3*(mu3^2+1-mu4)
        D = B^2 - A*(e3^2*mu3 + 2*e1*e3 - e4)
        
        if (D>=0){ 
          a2=(B+sqrt(D))/A
          a0 = e1 - a2
          a1 = e3 - a2*mu3
          tao2 = e2 - a0^2 - a1^2 - a2^2*mu4 - 2*a0*a2 - 2*a1*a2*mu3
          if (tao2 <= 0) return(NA)
          eta = matrix(rnorm(n*sim_num),nrow=n,ncol=sim_num)
          fgrs = (a0 + a1*y + a2*y^2) %*% t(rep(1,sim_num)) + sqrt(tao2)*eta
        }else{
          a0 = e1
          a1 = e3
          tao2 = e2- a0^2 - a1^2
          if (tao2 <=0) return(NA)
          eta = matrix(rnorm(n*sim_num),nrow=n,ncol=sim_num)
          fgrs = (a0 + a1*y) %*% t(rep(1,sim_num)) + sqrt(tao2)*eta
        }
        # return(fgrs)
        individual_coefficients_fgrs  =  lapply( as.data.frame(fgrs),
                                                 function(fg) {
                                                   .single_gxe( y,
                                                                as.matrix( fg ),
                                                                c( ao_het, 0.1, 0 ) )
                                                 } )
        
        individual_coefficients_fgrs =  do.call( cbind, individual_coefficients_fgrs )
        rownames( individual_coefficients_fgrs )  =  coef_names
        coefficients_fgrs    =  apply( individual_coefficients_fgrs, 1, mean )
        se_coefficients_fgrs =  apply( individual_coefficients_fgrs, 1, sd   )
        p_coefficients_fgrs  =  2 * pnorm( -abs( coefficients_fgrs / se_coefficients_fgrs ) )
        t_real_fgrs          =  ( coefficients - coefficients_fgrs ) / sqrt( se_coefficients^2 + se_coefficients_fgrs^2 )
        # results  =  list()
        
        return(c(list(es = c(e1=e1,e2=e2,e3=e3,e4=e4),
               fake_grs = list(
                coefficients = setNames( coefficients_fgrs, coef_names ),
                se = setNames( se_coefficients_fgrs, coef_names ),
                p  = setNames( p_coefficients_fgrs,  coef_names ),
                individual_estimates = individual_coefficients_fgrs
              ),
              # fake_phenotype = results,
              t_real_fgrs = setNames( t_real_fgrs,       coef_names ),
              fgrs = fgrs
        ) ) )
      }
      stopCluster(cl)
      
      results[[nn]]= list(
        real_data = list(
          coefficients = setNames( coefficients, coef_names ),
          se = setNames( se_coefficients, coef_names ),
          p  = setNames( p_coefficients,  coef_names ),
          individual_estimates = individual_coefficients
        ),
        fake_grs = fGRS
        # fake_phenotype = results,
      )
  }

save(results,file='./data/gamma_fG_2.Rdata')

############ plots #############
##for gamma_fG
results[[1]][["fake_grs"]] = results[[1]][["fake_grs"]][-c(1,12,111)]
##remove NA

allvar = list()
for (idx in 1:2){
  es = lapply(results[[idx]][["fake_grs"]],function(x) x[["es"]])
  es = do.call(rbind, es)
  pars = lapply(results[[idx]][["fake_grs"]],function(x) x[["fake_grs"]][["coefficients"]])
  pars = do.call(rbind, pars)
  
  var = outer(es[,3:4],c(1,2,3,4),function(x,y) x^y)
  size = dim(var)
  dim(var) = c(size[1],size[2]*size[3])
  colnames(var) = c('e3.1','e4.1','e3.2','e4.2','e3.3','e4.3','e3.4','e4.4')
  int = expand.grid(c('e3.1','e3.2','e3.3','e3.4'),c('e4.1','e4.2','e4.3','e4.4'))
  ##??? t(int)
  var_int = apply(t(int),2,function(x) var[,x[1]]*var[,x[2]]) ##???
  colnames(var_int) = paste(int[,1],int[,2],sep='_')
  var = cbind(gamma=pars[,4],var,var_int)
  var = as.data.frame(var)
  allvar[[idx]] = var
}

null = allvar[[1]]
inter = allvar[[2]]
# null = as.data.frame(apply(null,2,scale))
# inter = as.data.frame(apply(inter,2,scale))

res1 = lm(gamma~.,data=null)
summary(res1) ## 4.1 4.3 related
# res.aic = step(res1,trace=F)
# summary(res.aic)
res1 = lm(gamma~.,data=null)
res.bic = step(res1,k=log(nrow(null)),trace = F)
summary(res.bic)
var = colnames(res.bic$model)
inter2 = inter[,var]
res2 = lm(gamma~.,data=inter2)
summary(res2)
res2.bic = step(res2,k=log(nrow(null)),trace = F)
summary(res2.bic)

# coef = sort(abs(res1$coefficients))
res2 = lm(gamma~e4.1+e3.2+e3.4,data=null)
summary(res2)

null[,2:3] = round(null[,2:3],2)
inter[,2:3] = round(inter[,2:3],2)
ne3 = round(median(null$e3.1),2)
ne4 = round(median(null$e4.1),2)
ie3 = round(median(inter$e3.1),2)
ie4 = round(median(inter$e4.1),2)

xidx = seq(2,11,2)

pdf('./results/gamma_fG.pdf',width=10,height = 8)
par(mfrow=c(2,2))
plot(data=null,gamma~e3.1,
     xlab=paste0('e3: setted E(fG*Y) = ',ne3,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='No interaction, gamma=0')
plot(data=null,gamma~e4.1,
     xlab=paste0('e3: setted E(fG^2*Y) = ',ne4,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='No interaction, gamma=0')
plot(data=inter,gamma~e3.1,
     xlab=paste0('e3: setted E(fG*Y) = ',ie3,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='Interaction, gamma=sqrt(0.05')
plot(data=inter,gamma~e4.1, #xaxt='n',
     xlab=paste0('e3: setted E(fG^2*Y) = ',ie4,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='Interaction, gamma=sqrt(0.05')
# axis(side=1,at=round(unique(inter$e4.1)[xidx],2))
par(mfrow=c(2,2))
boxplot(data=null,gamma~e3.1,
     xlab=paste0('e3: setted E(fG*Y) = ',ne3,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='No interaction, gamma=0')
boxplot(data=null,gamma~e4.1,
     xlab=paste0('e3: setted E(fG^2*Y) = ',ne4,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='No interaction, gamma=0')
boxplot(data=inter,gamma~e3.1,
     xlab=paste0('e3: setted E(fG*Y) = ',ie3,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='Interaction, gamma=sqrt(0.05')
boxplot(data=inter,gamma~e4.1, #xaxt='n',
     xlab=paste0('e3: setted E(fG^2*Y) = ',ie4,'(real E) +- 0.25'),
     ylab='fake gamma',
     main='Interaction, gamma=sqrt(0.05')
dev.off()

############## try lm of diff gamma with grs ############
load('./data/fig4.Rdata')

## allcomb
alln    = seq(1e4,1e5,2e4) # sample size
allc1   = sqrt(seq(0.002,0.02,0.004)) # gamma
allcomb = expand.grid(alln,allc1)
colnames(allcomb) = c('sample','gamma')

idx = which(allcomb$sample==1e4)
comb1 = allcomb[idx,]
data = lapply(sim4,function(x) x[idx])

