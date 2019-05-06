
# purity calculation functions

# library
library(lme4)
library(mixtools)
library(MASS)

# main function, just run this 
purity_function = function(A, varname = '', times = '', 
                           trt = '',
                           trtlevel = '',
                           subj = '',
                           outcome = '',
                           start = 0, data = dat){
  p = length(A)
  alpha_est = matrix(A,p,1)
  dat_est = data
  if(sum(varname == '') == length(varname)){
    w_est = dat[(start):(start + p -1),] %*% alpha_est
  }
  if(sum(varname != '') == length(varname)){
    w_est = as.matrix(dat_est[,varname]) %*% alpha_est
  }
  dat_est$w = w_est
  if(times != ''){
    dat_est$tt = dat_est[,times]
  }
  if(subj != ''){
    dat_est$subj = dat_est[,subj]
  }
  if(trt != ''){
    dat_est$trt = dat_est[,trt]
    if(sum(trtlevel == '')==0){
      dat_est[dat_est$trt == trtlevel[1],]$trt = 'pbo'
      dat_est[dat_est$trt == trtlevel[2],]$trt = 'drg'
    }
  }
  dat_est$outcome = dat_est[,outcome]
  dat_pbo_est = dat_est[dat_est$trt == 'pbo', ]
  dat_drg_est = dat_est[dat_est$trt == 'drg', ]
  
  fit_drg_est = lmer(outcome ~ tt + I(tt^2) + w + w * tt +
                       w * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_drg_est, REML = FALSE)
  fit_drg_est
  fit_pbo_est = lmer(outcome ~ tt + I(tt^2) + w + w * tt +
                       w * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_pbo_est, REML = FALSE)
  fit_pbo_est
  
  beta1 = as.matrix(fixef(fit_drg_est))[2:3]
  gamma1 = as.matrix(fixef(fit_drg_est))[5:6] # true estimate = -1.9812346 -0.9895642
  D1 = as.matrix(VarCorr(fit_drg_est)$subj)[2:3, 2:3] 
  
  beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
  gamma2 = as.matrix(fixef(fit_pbo_est))[5:6] # 1.994051 1.003545
  D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[2:3, 2:3]
  
  b = purity_calculation(dat_est, beta1, beta2, gamma1, gamma2, D1, D2, p=2)
  return(list(purity = b$purity, Mu1 = b$Mu1, Mu2 = b$Mu2, beta1 = beta1, beta2 = beta2, 
              gamma1 = gamma1, gamma2 = gamma2, D1 = D1, D2 = D2, data = dat_est))
}

# calculate purity for each subject. 
purity_calculation = function(dat, beta1, beta2, gamma1, gamma2, D1, D2, p=2){
  unique_dat = unique(dat[,c('subj','w')])
  purity = c()
  Mu1 = c(); Mu2 = c()
  for(i in 1:dim(unique_dat)[1]){
    #if(i %% 10 ==0) print(i)
    mu1 = beta1 + gamma1 * unique_dat$w[i] # mu1 = beta_drg + gamma_drg * unique_dat$w[i]
    mu2 = beta2 + gamma2 * unique_dat$w[i] # mu2 = beta_pbo + gamma_pbo * unique_dat$w[i]
    Mu1 = rbind(Mu1, mu1); Mu2 = rbind(Mu2, mu2)
    res = mean(monta_carlo_pdf(Xstart, mu1, D1, mu2, D2))
    purity = c(purity, res)
  }
  return(list(purity = purity, Mu1 = Mu1, Mu2 = Mu2))
}

# Monte Carlo simulation
Xstart = mvrnorm(10000, c(0,0), diag(c(1,1)))
monta_carlo_pdf = function(Xstart, mu1, D1, mu2, D2){

  d1 = eigen(D1)
  d1 = d1$vectors %*% diag(sqrt(d1$values))
  
  d2 = eigen(D2)
  d2 = d2$vectors %*% diag(sqrt(d2$values))
  
  points1 = Xstart %*% t(d1) 
  points1[,1] = points1[,1] + mu1[1]; points1[,2] = points1[,2] + mu1[2]  
  points2 = Xstart %*% t(d2) 
  points2[,1] = points2[,1] + mu2[1]; points2[,2] = points2[,2] + mu1[2]  
  
  f1 = dmvnorm(points1, mu1, D1)
  f2 = dmvnorm(points2, mu2, D2)
  
  f3 = (f1 + f2) # deal with the 0s, f1 + f2 cannot be 0
  f3 = ifelse(f3 == 0, 1e-4,f3)
  
  purity = c((f1 - f2)^2 / f3) 
  return(purity)
}
