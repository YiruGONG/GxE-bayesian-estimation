OLS_1 = function(y,x,b0=TRUE){
  n = length(y)
  
  one = rep(1,n)
  xba = colMeans(x)
  XX  = x - one %*% t(xba)
  par1 = solve(t(XX)%*%XX) %*% t(XX) %*% y
  if (b0){
    yba = mean(y)
    par0 = yba - t(xba) %*% par1
    return(c(par0,par1))
  } else {
    return(par1)}
}

OLS_2 = function(y,x,par1){
  yr2 = (y - x%*%par1)^2
  x1  = x[,1]
  x2  = x[,2]
  
  lx1x1 = sum( (x1-mean(x1))^2 )
  lx2x2 = sum( (x2-mean(x2))^2 )
  lyy   = sum( (yr2-mean(yr2))^2 )
  lx1y  = sum( (x1-mean(x1))*(yr2-mean(yr2)) )
  lx2y  = sum( (x2-mean(x2))*(yr2-mean(yr2)) )
  lx1x2 = sum( (x1-mean(x1))*(x2-mean(x2)) )
  
  gamma2 = (lx2y*lx1x1 - lx1y*lx1x2) / (lx1x1*lx2x2 - lx1x2^2)
  if (gamma2<=0) {
    gamma=0
    func = function(par,y){ (mean(y)-par[1]^2-par[2])^2 }
    par2 = optim(c(0,0),func,y=yr2)$par
    beta = par2[1]
    sig2  = par2[2]
  } else {
    gamma=sqrt(gamma2)
    beta  = (lx1y - gamma^2*lx1x2) / (2*gamma*lx1x1)
    sig2   = mean(yr2) - 2*beta*gamma*mean(x1) - gamma^2*mean(x2) - beta^2
  }
  coeff = c(par1,beta,gamma,sig2)
  coeff
}
