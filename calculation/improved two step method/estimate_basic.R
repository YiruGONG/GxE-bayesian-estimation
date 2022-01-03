source('./calculation/two step/step1_fit.R')
source('./calculation/two step/step2_fit.R')
source('./calculation/two step/OLS.R')

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
    coef_names = c('alpha1', 'alpha2', 'beta', 'gamma','sigma' )
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
    coef_names = c('alpha1', 'alpha2', 'beta', 'gamma','sigma2' )
  } else if (method=='OLS.copar'){
    x1 = cbind(grs,(grs^2-1))
    par1 = lm(y~x1)$coefficients[-1]
    x2 = cbind(grs,grs^2)
    yr2 = (y - x1%*%par1)^2
    copar = lm(yr2~x2)$coefficients
    coeff = c(par1,copar)
    coef_names = c('alpha1', 'alpha2', 'b0', 'b1','b2' )
    # return(setNames( coeff, coef_names ))
  }
  return(setNames( coeff, coef_names )) 
}