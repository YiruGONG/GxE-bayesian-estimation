library(MASS) # simulate from  Multivariate Normal Distribution
library(LearnBayes) # generate inverse gamma sample
library(Directional) # vmf
library(RandomFieldsUtils) # solve linear equation for SPD matrix

# part 1: constraint prior
x_input = function(G, E)
{
  # G: n*1
  # E: n*1 missing values
  n = length(G)
  E = as.vector(E)
  p = 4
  x = matrix(0, ncol=p, nrow=n)
  x[,1] = G
  x[,2] = G*G-1
  x[,3] = E
  x[,4] = G*E
  return(x)
}
gram1 = function(x, beta, tau, l)
{
  # calculate the covariance matrix
  # x: n*4
  # beta: 4*1
  # tau,l: parameter
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
g_yebv1 = function(y, x, beta, var, tau, l)
{
  # y:n*1 x:n*p beta:p*1 var:number
  # tau,l:parameter
  n = length(y)
  p = length(beta)
  k = gram1(x, beta, tau, l)
  m = k + diag(rep(var, n))
  solution = solvex(m,k)     ##Cholesky decomposition
  miu_g = solution%*%y
  sigma_g = var*solution
  miu_g = as.vector(miu_g)
  g_new = mvrnorm(n=1, miu_g, sigma_g)
  return(g_new)
}
v_yegb1 = function(y, g, a=0.5, b=0.5)
{
  # y:n*1 g:n*1
  # a,b: hyperparameter in inv_gamma dist
  n = length(y)
  shape = a + n/2
  scale = b + 0.5*sum((y-g)^2)
  v_new = rigamma(n=1, shape, scale)
  return(v_new)
}
b_yegv1 = function(y, x, beta, var, tau, l, lamda=2000, var_b=10^(-3), sigma_b=1, proposal)
{
  # x : n*p beta:p*1 g:n*1 
  # p=4
  # lamda : vmf para
  # var_b : unif_proj/t_proj para
  # proposal = 'vmf','unif_proj','t_proj'
  # sigma_b: variance prior for beta 
  n = nrow(x)
  p = length(beta)

  if(proposal=='vmf')
  {
    beta_new = rvmf(n=1, mu=beta, k=lamda)
    beta_new = as.vector(beta_new)
  }
  else if (proposal=='unif_proj') {
    beta_new = mvrnorm(n=1, beta, diag(rep(var_b,p)))
    beta_new = as.vector(beta_new)
    if(sum(abs(beta_new)>1)>0){
      return(list(beta,0))
    }
  }
  else{
    beta_new = mvrnorm(n=1, beta, diag(rep(var_b,p)))
    beta_new = as.vector(beta_new)    
  }

  k = gram1(x, beta, tau, l)
  k_new = gram1(x, beta_new, tau, l)

  k_var = k + diag(rep(var,n))
  k_new_var = k_new + diag(rep(var,n))

  r = chol(k_var) # cholesky dcomposition
  r_new = chol(k_new_var) # cholesky dcomposition

  log_det_k_var = 2*sum(log(diag(r)))
  log_det_k_new_var = 2*sum(log(diag(r_new)))

  inv_l_y  = backsolve(r, y, transpose=TRUE)
  inv_l_new_y = backsolve(r_new, y, transpose=TRUE)

  log_ratio = 0.5*(log_det_k_var - log_det_k_new_var + crossprod(inv_l_y) - crossprod(inv_l_new_y) )[1,1]
  # if(proposal=='norm_proj')
  # {
  #   log_ratio = log_ratio + 0.5/sigma_b*(beta%*%beta - beta_new%*%beta_new)[1,1]
  # }
  if(proposal=='t_proj')
  {
    inv_sigma_b = diag(1/sigma_b)
    log_ratio = log_ratio + 0.5*(beta%*%inv_sigma_b%*%beta - beta_new%*%inv_sigma_b%*%beta_new)[1,1]
  }

  u = runif(1)
  if(log(u) <log_ratio){
    return(list(beta_new/sqrt(sum((beta_new)^2)), 1))
  }else{
    return(list(beta, 0))
  }
}
tau_yebgv1 = function(y, x, beta, tau, l,var, var_tau, a=0.5, b=0.5)
{
  n = nrow(x)
  p = length(beta)
  
  log_tau_new = rnorm(n=1, mean=log(tau), sd=sqrt(var_tau))
  tau_new = exp(log_tau_new)

  k = gram1(x, beta, tau, l)
  k_new = gram1(x, beta, tau_new, l)

  k_var = k + diag(rep(var,n))
  k_new_var = k_new + diag(rep(var,n))

  r = chol(k_var) # cholesky dcomposition
  r_new = chol(k_new_var) # cholesky dcomposition

  log_det_k_var = 2*sum(log(diag(r)))
  log_det_k_new_var = 2*sum(log(diag(r_new)))

  inv_l_y  = backsolve(r, y, transpose=TRUE)
  inv_l_new_y = backsolve(r_new, y, transpose=TRUE)

  log_ratio = 0.5*(log_det_k_var - log_det_k_new_var + crossprod(inv_l_y) - crossprod(inv_l_new_y) )[1,1]
  log_ratio = log_ratio-(a+1)*log(tau_new/tau) -b*(1/tau_new - 1/tau)
  log_ratio = log_ratio + log(tau_new) - log(tau)      ###???

  u = runif(1)
  if(log(u) <log_ratio){
    return(list(tau_new, 1))
  }else{
    return(list(tau, 0))
  }
}
l_yebgv1 = function(y, x, beta, tau, l, var, var_l, a=0.5, b=0.5)
{
  n = nrow(x)
  p = length(beta)
  
  log_l_new = rnorm(n=1, mean=log(l), sd=sqrt(var_l))
  l_new = exp(log_l_new)

  k = gram1(x, beta, tau, l)
  k_new = gram1(x, beta, tau, l_new)

  k_var = k + diag(rep(var,n))
  k_new_var = k_new + diag(rep(var,n))

  r = chol(k_var) # cholesky dcomposition
  r_new = chol(k_new_var) # cholesky dcomposition

  log_det_k_var = 2*sum(log(diag(r)))
  log_det_k_new_var = 2*sum(log(diag(r_new)))

  inv_l_y  = backsolve(r, y, transpose=TRUE)
  inv_l_new_y = backsolve(r_new, y, transpose=TRUE)

  log_ratio = 0.5*(log_det_k_var - log_det_k_new_var + crossprod(inv_l_y) - crossprod(inv_l_new_y) )[1,1]
  log_ratio = log_ratio-(a+1)*log(l_new/l) -b*(1/l_new - 1/l)
  log_ratio = log_ratio + log(l_new) - log(l)

  u = runif(1)
  if(log(u) <log_ratio){
    return(list(l_new, 1))
  }else{
    return(list(l, 0))
  }
}
sigma_b_yebgv1 = function(beta, a, b)
{
  p = length(beta)
  sigma_b_new = rep(0,p)
  for(i in 1:p)
  {
    shape = a + 1/2
    scale = b + 0.5*(beta[i])^2
    sigma_b_new[i] = rigamma(n=1, shape, scale)
  }
  return(sigma_b_new)
}
E_ybgv1 = function(y, G, E, beta, tau, l, var, var_E)
{
  n = length(y)
  p = length(beta)
  
  E = as.vector(E)
  # E_new = mvrnorm(n=1, E, diag(rep(var_E, n)))
  E_new = rnorm(n,mean(E),var_E)
  E_new = as.vector(E_new)

  x = x_input(G, E)
  x_new = x_input(G, E_new)

  k = gram1(x, beta, tau, l)
  k_new = gram1(x_new, beta, tau, l)

  k_var = k + diag(rep(var,n))
  k_new_var = k_new + diag(rep(var,n))

  r = chol(k_var) # cholesky dcomposition
  r_new = chol(k_new_var) # cholesky dcomposition

  log_det_k_var = 2*sum(log(diag(r)))
  log_det_k_new_var = 2*sum(log(diag(r_new)))

  inv_l_y  = backsolve(r, y, transpose=TRUE)
  inv_l_new_y = backsolve(r_new, y, transpose=TRUE)

  log_ratio = 0.5*(log_det_k_var - log_det_k_new_var + crossprod(inv_l_y) - crossprod(inv_l_new_y) )[1,1]
  log_ratio = log_ratio + 0.5*(sum(E^2) - sum(E_new^2))

  u = runif(1)
  if(log(u) <log_ratio){
    return(list(E_new, 1))
  }else{
    return(list(E, 0))
  }
}
mh_gib1 = function( y, G, E0, g0, beta0, var0, tau0, l0, lamda=2000, var_b=10^(-8), sigma_b0 = 0.5, a=0.5, b=0.5, var_E, var_tau, var_l, ite, proposal)
{
  # proposal = 'vmf','t_proj','unif_proj', #  'normal_proj'
  n = length(g0)
  p = length(beta0)

  E = matrix(0, ncol=n, nrow=(ite+1))
  g = matrix(0, ncol=n, nrow=(ite+1))
  beta = matrix(0, ncol=p, nrow=(ite+1))
  var = rep(0, ite+1)
  tau = rep(0, ite+1)
  l = rep(0, ite+1)
  
  E[1,] = E0
  g[1,] = g0
  beta[1,] = beta0
  var[1] = var0
  tau[1] = tau0
  l[1] = l0
  if(proposal=='t_proj')
  {
    sigma_b = matrix(0, ncol=p, nrow=(ite+1))
    sigma_b[1,] = sigma_b0
  }

  accept_E = 0
  accept_beta = 0
  accept_tau = 0
  accept_l = 0

  x = x_input(G, E0)
  for(i in 1:ite)
  {
    # if(proposal=='normal_proj')
    # {
    #   beta_outp = b_ygv1(x=x, beta=beta[i,], y, var=var[i], tau=tau[i], l=l[i], var_b=var_b, sigma_b=sigma_b0, proposal=proposal)
    # }
    # else
    if (proposal=='t_proj') {
      beta_outp = b_yegv1(y, x, beta[i,], var[i], tau[i], l[i], var_b=var_b, sigma_b=sigma_b[i,] ,proposal=proposal)
    }
    else {
      beta_outp = b_yegv1(y, x, beta[i,], var[i], tau[i], l[i], lamda=lamda, var_b=var_b, proposal=proposal)
    }
    
    beta[i+1,] = beta_outp[[1]]
    accept_beta = accept_beta + beta_outp[[2]]

    tau_outp = tau_yebgv1(y, x, beta[i+1,], tau[i], l[i], var[i], var_tau)
    tau[i+1] = tau_outp[[1]]
    accept_tau = accept_tau + tau_outp[[2]]

    l_outp = l_yebgv1(y, x, beta[i+1,], tau[i+1], l[i], var[i], var_l)
    l[i+1] = l_outp[[1]]
    accept_l = accept_l + l_outp[[2]] 
    
    E_outp = E_ybgv1(y, G, E[i,], beta[i+1,], tau[i+1], l[i+1], var[i], var_E)
    E[i+1,] = E_outp[[1]]
    accept_E = accept_E + E_outp[[2]]

    if(E_outp[[2]] == 1)
    {
      x = x_input(G, E[i+1,])
    }

    g[i+1,] = g_yebv1(y, x, beta[i+1,], var[i], tau[i+1], l[i+1])

    var[i+1] = v_yegb1(y, g[i+1,])

    if(proposal=='t_proj')
    {
      sigma_b[i+1,] = sigma_b_yebgv1(beta[i+1,], a, b)
    }
  }

  result = list(E, beta, g, var, tau, l, c(accept_beta,accept_E,accept_tau,accept_l)/ite)
  return(result)
}

#-------------------------------------------------------------------------------------------------------

# part 2: Unconstraint prior

# gram2 = function(beta, x, tau)
# {
#   # calculate the covariance matrix
#   # x:n*p, beta:p*1
#   # tau:parameter
#   covar = function(x1, x2, tau)
#   {
#     # GP covariance function
#     a = tau*exp(-(x1-x2)^2)
#     return(a)
#   }
#   n = nrow(x)
#   cov_m = matrix(0,n,n)
#   for (i in 1:n)
#   {
#     for (j in i:n)
#     {
#       if(i==j){
#         cov_m[i,i] = tau
#       } else{
#         xi = x[i,]%*%beta
#         xi = xi[1,1]
#         xj = x[j,]%*%beta
#         xj = xj[1,1]        
#         a = covar(xi, xj, tau)
#         cov_m[i,j] = a
#         cov_m[j,i] = a
#       }
#     }
#   }
#   return(cov_m)
# }
# g_ybv2 = function(y, x, beta, var, tau)
# {
#   # y:n*1 x:n*p beta:p*1 var:number
#   # tau:parameter
#   n = length(y)
#   p = length(beta)
#   k = gram2(beta=beta, x=x, tau=tau)
#   m = k + diag(rep(var, n))
#   solution = solvex(m,k)
#   miu_g = solution%*%y
#   sigma_g = var*solution
#   miu_g = as.vector(miu_g)
#   g_new = mvrnorm(n=1, miu_g, sigma_g)
#   return(g_new)
# }
# v_ygb2 = function(y, g, a=0.5, b=0.5)
# {
#   # y:n*1 g:n*1
#   # a,b: hyperparameter in inv_gamma dist
#   n = length(y)
#   shape = a + n/2
#   scale = b + 0.5*sum((y-g)^2)
#   v_new = rigamma(n=1, shape, scale)
#   return(v_new)
# }
# sigma_b_ybgv2 = function(beta, a, b)
# {
#   p = length(beta)
#   sigma_b_new = rep(0,p)
#   for(i in 1:p)
#   {
#     shape = a + 1/2
#     scale = b + 0.5*(beta[i])^2
#     sigma_b_new[i] = rigamma(n=1, shape, scale)
#   } 
#   return(sigma_b_new)
# }
# b_ygv2 = function(x, beta, y, var, tau, sigma_b, var_b)
# {
#   # x:n*p beta:p*1
#   # var_b:normal random walk para
#   # sigma_b: prior variance for beta
#   n = nrow(x)
#   p = length(beta)

#   beta_new = mvrnorm(n=1, beta, diag(rep(var_b,p)))
#   beta_new = as.vector(beta_new)

#   k = gram2(beta=beta, x=x, tau=tau)
#   k_new = gram2(beta=beta_new, x=x, tau=tau)

#   k_var = k + diag(rep(var,n))
#   k_new_var = k_new + diag(rep(var,n))

#   r = chol(k_var) # cholesky dcomposition
#   r_new = chol(k_new_var) # cholesky dcomposition

#   log_det_k_var = 2*sum(log(diag(r)))
#   log_det_k_new_var = 2*sum(log(diag(r_new)))

#   inv_l_y  = backsolve(r, y, transpose=TRUE)
#   inv_l_new_y = backsolve(r_new, y, transpose=TRUE)

#   log_ratio = 0.5*(log_det_k_var - log_det_k_new_var + crossprod(inv_l_y) - crossprod(inv_l_new_y) )[1,1]
#   inv_sigma_b = diag(1/sigma_b)
#   log_ratio = log_ratio + 0.5*(beta%*%inv_sigma_b%*%beta - beta_new%*%inv_sigma_b%*%beta_new)[1,1]

#   u = runif(1)
#   if(log(u) <log_ratio){
#     return(list(beta_new, 1))
#   }else{
#     return(list(beta, 0))
#   }
# }
# tau_ybgv2 = function(beta, y, x, tau,var, var_tau, a=0.5, b=0.5)
# {
#   n = nrow(x)
#   p = length(beta)
  
#   log_tau_new = rnorm(n=1, mean=log(tau), sd=sqrt(var_tau))
#   tau_new = exp(log_tau_new)

#   k = gram2(beta=beta, x=x, tau=tau)
#   k_new = gram2(beta=beta, x=x, tau=tau_new)

#   k_var = k + diag(rep(var,n))
#   k_new_var = k_new + diag(rep(var,n))

#   r = chol(k_var) # cholesky dcomposition
#   r_new = chol(k_new_var) # cholesky dcomposition

#   log_det_k_var = 2*sum(log(diag(r)))
#   log_det_k_new_var = 2*sum(log(diag(r_new)))

#   inv_l_y  = backsolve(r, y, transpose=TRUE)
#   inv_l_new_y = backsolve(r_new, y, transpose=TRUE)

#   log_ratio = 0.5*(log_det_k_var - log_det_k_new_var + crossprod(inv_l_y) - crossprod(inv_l_new_y) )[1,1]
#   log_ratio = log_ratio-(a+1)*log(tau_new/tau) -b*(1/tau_new - 1/tau)
#   log_ratio = log_ratio + log(tau_new) - log(tau)
  
#   u = runif(1)
#   if(log(u) <log_ratio){
#     return(list(tau_new, 1))
#   }else{
#     return(list(tau, 0))
#   }
# }
# mh_gib2 = function(g0, beta0, var0, y, x, tau0, sigma_b0, a, b, var_b, var_tau, ite)
# {
#   #
#   n = length(g0)
#   p = length(beta0)

#   g = matrix(0, ncol=n, nrow=(ite+1))
#   beta = matrix(0, ncol=p, nrow=(ite+1))
#   var = rep(0, ite+1)
#   tau = rep(0, ite+1)
#   sigma_b = matrix(0, ncol=p, nrow=(ite+1))
  
#   accept_beta = 0
#   accept_tau = 0
#   log_u = rep(0,ite)
#   log_ratio = rep(0,ite)

#   g[1,] = g0
#   beta[1,] = beta0
#   var[1] = var0
#   tau[1] = tau0
#   sigma_b[1,] = sigma_b0
#   # gibbs sampling
#   for(i in 1:ite)
#   {
#     beta_outp = b_ygv2(x, beta[i,], y, var[i], tau[i], sigma_b[i,], var_b=var_b)
#     beta[i+1,] = beta_outp[[1]]
#     accept_beta = accept_beta + beta_outp[[2]]

#     g[i+1,] = g_ybv2(y=y, x=x, beta=beta[i+1,], var=var[i], tau=tau[i])
    
#     var[i+1] = v_ygb2(y=y, g=g[i+1,])

#     sigma_b[i+1,] = sigma_b_ybgv2(beta[i+1,],a=a,b=b)

#     tau_outp = tau_ybgv2(beta[i+1,], y, x, tau[i], var[i+1], var_tau)
#     tau[i+1] = tau_outp[[1]]
#     accept_tau = accept_tau + tau_outp[[2]]
#   }

#   result = list(g,var,beta,tau,sigma_b,c(accept_beta,accept_tau)/ite)
#   return(result)
# }


