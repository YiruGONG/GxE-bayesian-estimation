fit_1 = function(par,y,grs){
  a = par[1]
  beta = par[2]
  gamma = par[3]
  sig2 = 1 - a^2 - beta^2 - gamma^2
  
  if (sig2 <= 0) {func = 0.5*sum(y^2)} else{
    miu = a*grs
    sigma2 = (beta+gamma*grs)^2 + sig2
    func = 0.5 * sum( (y-miu)^2 / sigma2 + log(sigma2) )
  }
  return(func)
}

fit_2 = function(par,y,grs){
  a1 = par[1]
  a2 = par[2]
  beta = par[3]
  gamma = par[4]
  sig2 = 1 - a1^2 - 2*a2^2 - beta^2 - gamma^2
  if (sig2 <= 0) {func = 0.5*sum(y^2)} else{
    miu = a1*grs + a2*(grs^2-1)
    sigma2 = (beta+gamma*grs)^2 + sig2
    func = 0.5 * sum( (y-miu)^2 / sigma2 + log(sigma2) )
  }
  return(func)
}

est_par = function(par,y,grs){
  par = optim(par,fit_1,y=y,grs=grs)$par
  optim(c(par[1],0,par[-1]),fit_2,y=y,grs=grs)$par
}

try_fG = function(y,grs,sim_num){
  mu1 = mean(y)
  mu2 = mean(y^2)
  mu3 = mean(y^3)
  mu4 = mean(y^4)
  mu5 = mean(y^5)
  b1 = mean(y*grs)
  b2 = mean(y*grs^2)
  
  A = mu3^3 + mu5 - 2*mu3*mu4
  B = - (b1*(mu4-mu3^2-1))
  Delta = B^2 - A * (b1^2*mu3 - b2)
  noi = t(rnorm(sim_num) %*% matrix(1,nrow=1, ncol=length(y)))
  
  if (Delta > 0){
    a2 = (B + sqrt(Delta))/A
    a0 = -a2 * mu2
    a1 = b1 - a2*mu3
    tao2 = 1 - 3*a2^2 + a1^2 + a2^2*mu4 + 2*a1*a2*mu3
    # tao2 = 1 - a0^2 - a1^2 - a2^2*mu4 - 2*a0*a2 - 2*a1*a2*mu3
    if (tao2 < 0) stop('tao2 < 0, stop on boundary')
    G0 = a0 + a1*y + a2*y^2
  } else{
    a0 = 0
    a1 = b1/miu2
    tao2 = 1 - a0^2 - a1^2*miu2
    if (tao2 < 0) stop('tao2 < 0, stop on boundary')
    G0 = a0 + a1*y
  }
  fG = G0 %*% matrix(1,nrow=1,ncol=sim_num) + sqrt(tao2)*noi
  return(fG)
}

find_best_a1 = function(y,grs,skws,krt,a1){
  all_a1 = seq(a1-0.25,a1+0.25,0.01)
  noi = rpearson( length(y), moments = c( 0, 1, skw, krt ) )
  z = grs %*% t(all_a1) + noi %*% t(sqrt(1-all_a1^2))
  ys = apply(z,2,function(z1){
    poly = outer(z1,1:max_polynomial,function(z,x) z^x)
    regpar = lm(y~poly)$coefficients
    y1 = scale(poly %*% regpar[2:8] + regpar[1])
    return(y1)
  })
  comp_a1 = apply(ys,2,function(y1,grs) lm(y1~grs)$coefficients[2],grs)
  na1 = which.min(abs(comp_a1-a1))
  # fa1 = comp_a1[na1]
  # fy = y[,na1]
  return(comp_a1[na1])
}


try_gxe = function(y,grs,sim_num, generate_fY = FALSE,
                   skws = seq(-3,3,0.2),
                   kurt_range = c(2:4),
                   max_polynomial = 7){
  g2 = grs^2
  a1 = lm(y~grs+g2)$coefficients[2]
  allpar = sapply(1:sim_num,function(x){
    idx = sample(1:length(y),length(y),replace = T)
    est_par(c(a1,0.1,0),y[idx],grs[idx])
  })
  parnames = c('a1','a2','beta','gamma')
  rownames(allpar) = parnames
  finalpar = rowMeans(allpar)
  
  ###generate fake GRS
  fG = try_fG(y,grs,sim_num)
  allpar_fG = sapply(as.data.frame(fG),function(fg) 
    est_par(c(a1,0.1,0),y,fg))
  rownames(allpar_fG) = parnames
  
  ###generate fake Y
  if (generate_fY) {
    momcomb = expand.grid(skw = skws, kurt_r = kurt_range)
    momcomb$krt = momcomb$skw^2 + momcomb$kurt_r
    a1 = as.numeric(cor(y,grs))
    origpar = optim(c(a1,0.1,0),fit_1,y=y,grs=grs)$par
    names(origpar) = c('alpha','beta','gamma')
    a1 = origpar[1]
    # allcomb = expand.grid(momidx = 1: nrow(momcomb),a1 = all_a1)
    res.comb = apply(momcomb,1,function(comb){
      skw = momcomb[comb[1],1]
      krt = momcomb[comb[1],3]
      best.a1 = find_best_a1(y,grs,skw,krt,a1)
      noi = matrix(rpearson( length(y)*sim_num, moments = c( 0, 1, skw, krt ) ),ncol=sim_num)
      zs = grs %*% matrix(best.a1,nrow=1,ncol=sim_num) + sqrt(1-best.a1^2)*noi
      ys = apply(zs,2,function(z1) {
        poly = outer(z1,1:max_polynomial,function(z,x) z^x)
        regpar = lm(y~poly)$coefficients
        y1 = scale(poly %*% regpar[2:8] + regpar[1])
        return(y1)
      }
      allpar_fY = apply(ys,2,function(y1) optim(origpar,fit_1,y=y1,grs=grs)$par))
      
      meanpar = rowMeans(allpar_fY)
      se = apply(allpar_fY,1,sd) / sqrt(sim_num)
      score = sum(abs(meanpar-origpar)/se)
      
      list(meanpar = meanpar,
           se = se,
           p = 2 * pnorm( -abs( coefficients / se ) )
           score = score,
           skwness = skw,
           kurtosis = krt,
           best.a1 = best.a1)
    })
    score = sapply(res.comb,function(x) x$score)
    results = res.comb[which.min(score)]
  }
  
  
}



