
# data generation for simulation
# can work for high dimensions

true_generation = function(alpha, p, n, ni, tt, X, 
                           beta_drg, gamma_drg, bi_sigma, sigma_drg,
                           beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo){
  # alpha
  set.seed(123)
  alpha = as.matrix(alpha,p,1)
  dat_drg = c()
  for(i in 1:n){
    drg_temp = NULL
    drg_temp$subj = rep(paste('drg',i,sep=''),ni)
    drg_temp$trt = rep('drg',ni)
    drg_temp = as.data.frame(drg_temp)
    baseline = as.matrix(rnorm(p,0,1),p,1)
    drg_temp = cbind(matrix(rep(baseline, each = ni),ni,p),drg_temp)
    colnames(drg_temp)[1:p] = paste('X',1:p, sep = '')
    w = rep(t(alpha) %*% baseline,ni)
    drg_temp$w = w
    drg_temp$tt = tt
    bi = mvrnorm(1, c(0,0,0), bi_sigma)
    yi = X%*%(beta_drg+bi+gamma_drg*w[1]) + sigma_drg*rnorm(ni,0,1)
    drg_temp$y = yi
    dat_drg = rbind(dat_drg, as.data.frame(drg_temp))
  }
  
  dat_pbo = c()
  for(i in 1:n){
    pbo_temp = NULL
    pbo_temp$subj = rep(paste('pbo',i,sep=''),ni)
    pbo_temp$trt = rep('pbo',ni)
    pbo_temp = as.data.frame(pbo_temp)
    baseline = as.matrix(rnorm(p),p,1)
    pbo_temp = cbind(matrix(rep(baseline, each = ni),ni,p),pbo_temp)
    colnames(pbo_temp)[1:p] = paste('X',1:p, sep = '')
    w = rep(t(alpha) %*% baseline,ni)
    pbo_temp$w = w
    pbo_temp$tt = tt
    bi = mvrnorm(1, c(0,0,0), (bi_sigma2))
    yi = X%*%(beta_pbo+bi+gamma_pbo*w[1]) + sigma_pbo*rnorm(ni,0,1)
    pbo_temp$y = yi
    dat_pbo = rbind(dat_pbo, as.data.frame(pbo_temp))
  }
  print('True data generated')
  return(list(dat_drg = dat_drg, dat_pbo = dat_pbo))
}