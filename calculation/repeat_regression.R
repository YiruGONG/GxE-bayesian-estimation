library(MASS) # simulate from  Multivariate Normal Distribution
library(LearnBayes) # generate inverse gamma sample
library(Directional) # vmf
library(RandomFieldsUtils) # solve linear equation for SPD matrix

# ##for test
# x = cbind(rnorm(100),runif(100))
# beta0 = c(0.1,-0.3)
# gf = function(a) a^2
# g0 = gf(x%*%beta0)
# sig0 = 0.5
# tau0 = 5
# l0 = 0.1
# y = g0 + sig0*rnorm(100)

core_k = function(x,beta,tau,l){
  # n = nrow(x)
  # k = matrix(nrow=n,ncol=n)
  # beta = as.vector(beta)
  # for (i in 1:n){
  #   for (j in 1:n){
  #     k[i,j]= tau*exp( - (x[i,]%*%beta - x[j,]%*%beta)^2 / l )[1,1]
  #   }
  # }
  # return(k)
  covar = function(x1, x2, tau, l)
  {
    # GP covariance function
    a = tau*exp(-(x1-x2)^2/l)
    return(a)
  }
  n = nrow(x)
  cov_m = matrix(0, n, n)
  for (i in 1:n)
  {
    for (j in i:n)
    {
      if(i==j){
        cov_m[i,i] = tau
      } else{
        xi = x[i,]%*%beta
        xi = xi[1,1]
        xj = x[j,]%*%beta
        xj = xj[1,1]        
        a = covar(xi, xj, tau, l)
        cov_m[i,j] = a
        cov_m[j,i] = a
      }
    }
  }
  return(cov_m)
}

beta_ite= function(y,x,beta,tau,l,g,sig,lambda=2000,method,sig_b,a,b){        ##vmf method
  ## method = vmf or unif or t_proj
  n = length(g)
  p = ncol(x)
  
  if (method=='vmf'){
    beta_new = rvmf(1,mu=beta,k=lambda)[1,] ###what if beta is single number
  } else if (method=='unif'){
    beta_new = runif(p,-1,1)
  } else if (method =='t_proj'){
    beta_new = mvrnorm(p,rep(0,4),diag(sig_b))
  }

  k = core_k(x,beta,tau,l)
  k_new = core_k(x,beta_new,tau,l)
  
  L = chol(k + sig*diag(nrow=n)) # cholesky dcomposition
  L_new =  chol(k_new + sig*diag(nrow=n)) # cholesky dcomposition
  
  log_det = 2*sum(log(diag(L)))
  log_det_new = 2*sum(log(diag(L_new)))
  
  inv_var = backsolve(L, y, transpose=TRUE)
  inv_var_new = backsolve(L_new, y, transpose=TRUE)
  
  log_ratio = 0.5*(log_det-log_det_new)+
    0.5*(crossprod(inv_var)-crossprod(inv_var_new))
  if (method=='t_proj'){
    inv_cov = diag(1/sig_b)
    log_ratio = log_ratio + 0.5*(t(beta)%*%inv_cov%*%beta - t(beta_new)%*%inv_cov%*%beta_new)
  }
  
  u = runif(1)
  if(log(u)<log_ratio) {
    return(list(beta_new/sqrt(sum((beta_new)^2)),1))
  }else { return(list(beta,0)) }
}

tau_ite = function(y,x,beta,tau,l,g,sig,var_tau,a=0.5,b=0.5){
  n = length(g)
  p = ncol(x)
  
  log_tau_new = rnorm(1,log(tau),sqrt(var_tau))
  tau_new = exp(log_tau_new)
  
  k = core_k(x,beta,tau,l)
  k_new = core_k(x,beta,tau_new,l)
  
  L = chol(k + sig*diag(nrow=n)) # cholesky dcomposition
  L_new =  chol(k_new + sig*diag(nrow=n)) # cholesky dcomposition
  
  # k_new_var = k_new + diag(rep(sig,n))
  # r_new = chol(k_new_var) # cholesky dcomposition
  
  log_det = 2*sum(log(diag(L)))
  log_det_new = 2*sum(log(diag(L_new)))
  
  inv_var = backsolve(L, y, transpose=TRUE)
  inv_var_new = backsolve(L_new, y, transpose=TRUE)
  
  log_ratio = 0.5*(log_det-log_det_new)+
    0.5*(crossprod(inv_var)-crossprod(inv_var_new))[1,1]+
    (a+1)*(log(tau)-log(tau_new)) + b*(1/tau-1/tau_new)+
    log(tau_new) - log(tau)
  
  u = runif(1)
  if(log(u)<log_ratio) {
    return(list(tau_new,1))
  }else { return(list(tau,0)) }
}

l_ite = function(y,x,beta,tau,l,g,sig,var_l,a,b){
  n = length(g)
  p = ncol(x)
  
  log_l_new = rnorm(1,log(l),sqrt(var_l))
  l_new = exp(log_l_new)
  
  k = core_k(x,beta,tau,l)
  k_new = core_k(x,beta,tau,l_new)
  
  L = chol(k + sig*diag(nrow=n)) # cholesky dcomposition
  L_new =  chol(k_new + sig*diag(nrow=n)) # cholesky dcomposition
  
  log_det = 2*sum(log(diag(L)))
  log_det_new = 2*sum(log(diag(L_new)))
  
  inv_var = backsolve(L, y, transpose=TRUE)
  inv_var_new = backsolve(L_new, y, transpose=TRUE)
  
  log_ratio = 0.5*(log_det-log_det_new)+
    0.5*(crossprod(inv_var)-crossprod(inv_var_new))+
    (a+1)*(log(l)-log(l_new)) + b*(1/l-1/l_new)+
    log(l_new) - log(l)
  
  u = runif(1)
  if(log(u)<log_ratio) {
    return(list(l_new,1))
  }else { return(list(l,0)) }
}

g_ite = function(y,x,beta,tau,l,g,sig){
  n = length(g)
  p = ncol(x)
  
  k = core_k(x,beta,tau,l)
  k_var = k + sig*diag(nrow=n)
  sol = solvex(k_var,k)
  cov = sig * sol
  miu = as.vector(sol %*% y)
  
  g_new = mvrnorm(1,miu,cov)
  return(g_new)
}

sig_ite = function(y,g,a=0.5,b=0.5){
  n = length(g)
  
  a1 = a+n/2
  # b1 = 0.5*(t(y-g)%*%(y-g)) + b
  b1 = b + 0.5*sum((y-g)^2)
  
  sig_new = rigamma(1,a1,b1)
  return(sig_new)
}

###general estimation
estimate_bayes = function(y,x,beta0,tau0,l0,g0,sig0,var_tau,var_l,a,b,iteration=2000,method='vmf'){
  g0 = as.vector(g0)
  n = length(g0)
  p = ncol(x)
  
  beta = matrix(nrow=iteration+1,ncol=p)
  g <- matrix(nrow=iteration+1,ncol=n)
  tau <- l <- sig <- rep(NA,iteration+1)
  
  beta[1,] <- beta0
  g[1,] <- g0
  tau[1] <- tau0
  l[1] <- l0
  sig[1] <- sig0
  
  acc_beta <- acc_tau <- acc_l <- 0
  
  for (i in 1:iteration){
    beta1 = beta_ite(y,x,beta[i,],tau[i],l[i],g[i,],sig[i],lambda=2000,method)
    beta[i+1,] = beta1[[1]]
    acc_beta = acc_beta + beta1[[2]]
    
    tau1 = tau_ite(y,x,beta[i+1,],tau[i],l[i],g[i,],sig[i],var_tau,a,b)
    tau[i+1] = tau1[[1]]
    acc_tau = acc_tau +tau1[[2]]
    
    l1 = l_ite(y,x,beta[i+1,],tau[i+1],l[i],g[i,],sig[i],var_l,a,b)
    l[i+1] = l1[[1]]
    acc_l = acc_l + l1[[2]]
    
    g[i+1,] = g_ite(y,x,beta[i+1,],tau[i+1],l[i+1],g[i,],sig[i])
    
    sig[i+1] = sig_ite(y,g[i+1,],a,b)
  }
  return(list(
    beta = beta,
    g = g,
    tau = tau,
    l = l,
    sigma2 = sig,
    accept_rate = c(beta=acc_beta, tau=acc_tau, l=acc_l)/iteration
  ))
}