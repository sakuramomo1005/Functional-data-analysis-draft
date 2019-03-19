## simulation with monte carlo
## date 2019-03-17


# set parameters:
# generate beta randomly
beta_drg = as.matrix(c(0,-5,-1),3,1)
beta_pbo = as.matrix(c(0,5,1),3,1)

# generate gamma randomly
gamma_drg=matrix(c(0,-2,-1),3,1)
gamma_pbo=matrix(c(0,2,1),3,1)

# bi
set.seed(123)
Aa = matrix(runif(3^2)*2-1, ncol=3) 
bi_sigma = t(Aa) %*% Aa
bi_sigma

# sigma
sigma_drg = 1
sigma_pbo = 1

### data generation
true_generation = function(alpha){
  # alpha
  set.seed(123)
  alpha = as.matrix(alpha,p,1)
  
  dat_drg = c()
  for(i in 1:n){
    drg_temp = NULL
    drg_temp$subj = rep(paste('drg',i,sep=''),ni)
    drg_temp$trt = rep('drg',ni)
    baseline = as.matrix(rnorm(p,0,1),p,1)
    x1 = baseline[1]; x2 = baseline[2]
    w = rep(t(alpha) %*% baseline,ni)
    drg_temp$x1 = rep(x1,ni); drg_temp$x2 = rep(x2,ni)
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
    baseline = as.matrix(rnorm(p),p,1)
    x1 = baseline[1]; x2 = baseline[2]
    w = rep(t(alpha) %*% baseline,ni)
    pbo_temp$x1 = rep(x1,ni); pbo_temp$x2 = rep(x2,ni)
    pbo_temp$w = w
    pbo_temp$tt = tt
    bi = mvrnorm(1, c(0,0,0), bi_sigma)
    yi = X%*%(beta_pbo+bi+gamma_pbo*w[1]) + sigma_pbo*rnorm(ni,0,1)
    pbo_temp$y = yi
    dat_pbo = rbind(dat_pbo, as.data.frame(pbo_temp))
  }
  
  return(list(dat_drg = dat_drg, dat_pbo = dat_pbo))
}

tt = as.matrix(0:6) # pt = the order of time points
ni = length(tt) # 7
X = cbind(matrix(1, length(tt), 1), tt, tt^2)

n = 2000 # number of subjects in each group  
p = 2 # number of baseline covariates
ni = 7 # number of time points

theta_true = 60/180 * pi # the true theta
data = true_generation(c(sin(theta_true),cos(theta_true)))
dat_drg = data$dat_drg
dat_pbo = data$dat_pbo
dat = rbind(dat_drg, dat_pbo)
head(dat)
dim(dat)

## check the estimates
fit_drg_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_drg, REML = FALSE)
fit_pbo_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_pbo, REML = FALSE)

beta1 = as.matrix(fixef(fit_drg_est))[2:3]
gamma1 = as.matrix(fixef(fit_drg_est))[5:6]
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[2:3, 2:3]

beta1; beta_drg; gamma1; gamma_drg; D1; bi_sigma;

beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
gamma2 = as.matrix(fixef(fit_pbo_est))[5:6]
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[2:3, 2:3]

beta2; beta_pbo; gamma2; gamma_pbo; D2; bi_sigma


# purity calculation
purity_calculation = function(dat, beta1, beta2, gamma1, gamma2, D1, D2, p=2){
  unique_dat = unique(dat[,c('subj','w')])
  purity = c()
  for(i in 1:dim(unique_dat)[1]){
    #if(i %% 10 ==0) print(i)
    mu1 = beta1 + gamma1 * unique_dat$w[i] # mu1 = beta_drg + gamma_drg * unique_dat$w[i]
    mu2 = beta2 + gamma2 * unique_dat$w[i] # mu2 = beta_pbo + gamma_pbo * unique_dat$w[i]
    res = monta_carlo_pdf(Xstart, mu1, D1, mu2, D2)/1000
    purity = c(purity, res)
  }
  return(purity)
}
purity_function = function(A, data = dat){
  alpha_est = matrix(A,2,1)
  dat_est = dat
  w_est = cbind(dat$x1,dat$x2) %*% alpha_est
  dat_est$w = w_est
  dat_pbo_est = dat_est[dat_est$trt == 'pbo', ]
  dat_drg_est = dat_est[dat_est$trt == 'drg', ]
  
  fit_drg_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                       w * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_drg_est, REML = FALSE)
  fit_drg_est; beta_drg; gamma_drg
  fit_pbo_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                       w * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_pbo_est, REML = FALSE)
  fit_pbo_est; beta_pbo; gamma_pbo
   
  beta1 = as.matrix(fixef(fit_drg_est))[2:3]
  gamma1 = as.matrix(fixef(fit_drg_est))[5:6] # true estimate = -1.9812346 -0.9895642
  D1 = as.matrix(VarCorr(fit_drg_est)$subj)[2:3, 2:3] 
  
  beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
  gamma2 = as.matrix(fixef(fit_pbo_est))[5:6] # 1.994051 1.003545
  D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[2:3, 2:3]
  
  b = purity_calculation(dat_est, beta1, beta2, gamma1, gamma2, D1, D2, p=2)
  return((b))
}
monta_carlo_pdf = function(Xstart, mu1, D1, mu2, D2){
  
  n_mc = dim(Xstart)[1]
  mu1 = matrix(rep(mu1, n_mc),n_mc,2,byrow = TRUE)
  mu2 = matrix(rep(mu2, n_mc),n_mc,2,byrow = TRUE)
  
  Q1 = diag((-1/2)*(Xstart-mu1)%*%solve(D1)%*%t(Xstart-mu1))
  Q2 = diag((-1/2)*(Xstart-mu2)%*%solve(D2)%*%t(Xstart-mu2))
  
  f1 = (1/(2*pi))*(1/sqrt(det(D1)))*exp(Q1)
  f2 = (1/(2*pi))*(1/sqrt(det(D2)))*exp(Q2)
  
  f3 = (f1 + f2) # deal with the 0s, f1 + f2 cannot be 0
  f3 = ifelse(f3 == 0, 1e-4,f3)
  
  purity = c((f1 - f2)^2 / f3) 
  return(purity)
}

# check purity value

# how sensitive the purity is?
# mu1 = c(3.523528, 0.2425273)
# mu2 = c(-3.493393, -0.2547145)
# D1
# D2
sum(monta_carlo_pdf(Xstart, mu1, D1, mu2, D2)) # 11.1564
sum(monta_carlo_pdf(Xstart, mu1+0.001, D1, mu2, D2)) # 11.19739
sum(monta_carlo_pdf(Xstart, mu1, D1, mu2, D2+0.001)) # 11.12985

# true purity
truep = sum(purity_calculation(dat, c(beta_drg)[2:3], c(beta_pbo)[2:3], 
                       c(gamma_drg)[2:3], c(gamma_pbo)[2:3], 
                       bi_sigma[2:3,2:3], bi_sigma[2:3,2:3], p=2))
fit_drg_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_drg, REML = FALSE)
fit_pbo_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_pbo, REML = FALSE)
beta1 = as.matrix(fixef(fit_drg_est))[2:3]
gamma1 = as.matrix(fixef(fit_drg_est))[5:6]
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[2:3, 2:3]
beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
gamma2 = as.matrix(fixef(fit_pbo_est))[5:6]
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[2:3, 2:3]

# any random purity
true_est_p = sum(purity_calculation(dat, beta1, beta2, gamma1, gamma2, D1, D2, p=2))

A = c(sin(theta_true),cos(theta_true))

unique_dat = unique(dat[,c('subj','w')])
mu_drg = c(); mu_pbo = c(); mu1_est = c(); mu2_est = c()
for(i in 1:dim(unique_dat)[1]){
  #if(i %% 10 ==0) print(i)
  mu1 = beta1 + gamma1 * unique_dat$w[i] # mu1 = beta_drg + gamma_drg * unique_dat$w[i]
  mu2 = beta2 + gamma2 * unique_dat$w[i] # mu2 = beta_pbo + gamma_pbo * unique_dat$w[i]
  mu1_est = c(mu1_est, mu1)
  mu2_est = c(mu2_est, mu2)
  mu_drg = c(mu_drg, beta_drg[2:3] + gamma_drg[2:3] * unique_dat$w[i])
  mu_pbo = c(mu_pbo, beta_pbo[2:3] + gamma_pbo[2:3] * unique_dat$w[i])
}
mus = data.frame(mu1_est = mu1_est, mu_drg = mu_drg,mu2_est = mu2_est, mu_pbo = mu_pbo)

mean(abs(mus$mu1_est - mus$mu2_est)) # 6.176789
mean(abs(mus$mu_drg - mus$mu_pbo)) # 6.172902

# how about a very different A?
A = c(1,1)
alpha_est = matrix(A,2,1)
dat_est = dat
w_est = cbind(dat$x1,dat$x2) %*% alpha_est
dat_est$w = w_est
dat_pbo_est = dat_est[dat_est$trt == 'pbo', ]
dat_drg_est = dat_est[dat_est$trt == 'drg', ]

fit_drg_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_drg_est, REML = FALSE)
fit_drg_est; beta_drg; gamma_drg
fit_pbo_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_pbo_est, REML = FALSE)
fit_pbo_est; beta_pbo; gamma_pbo

beta1 = as.matrix(fixef(fit_drg_est))[2:3]
gamma1 = as.matrix(fixef(fit_drg_est))[5:6] # true estimate = -1.9812346 -0.9895642
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[2:3, 2:3] 

beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
gamma2 = as.matrix(fixef(fit_pbo_est))[5:6] # 1.994051 1.003545
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[2:3, 2:3]

unique_dat = unique(dat_est[,c('subj','w')])
mu1_est = c(); mu2_est = c()
for(i in 1:dim(unique_dat)[1]){
  #if(i %% 10 ==0) print(i)
  mu1 = beta1 + gamma1 * unique_dat$w[i] # mu1 = beta_drg + gamma_drg * unique_dat$w[i]
  mu2 = beta2 + gamma2 * unique_dat$w[i] # mu2 = beta_pbo + gamma_pbo * unique_dat$w[i]
  mu1_est = c(mu1_est, mu1)
  mu2_est = c(mu2_est, mu2)
}
mus = data.frame(mu1_est = mu1_est, mu2_est = mu2_est)

mean(abs(mus$mu1_est - mus$mu2_est)) # 6.156634



A = c(20,0)
alpha_est = matrix(A,2,1)
dat_est = dat
w_est = cbind(dat$x1,dat$x2) %*% alpha_est
dat_est$w = w_est
dat_pbo_est = dat_est[dat_est$trt == 'pbo', ]
dat_drg_est = dat_est[dat_est$trt == 'drg', ]

fit_drg_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_drg_est, REML = FALSE)
fit_drg_est; beta_drg; gamma_drg
fit_pbo_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                     w * I(tt^2) + (tt+I(tt^2)|subj),
                   data = dat_pbo_est, REML = FALSE)
fit_pbo_est; beta_pbo; gamma_pbo

beta1 = as.matrix(fixef(fit_drg_est))[2:3]
gamma1 = as.matrix(fixef(fit_drg_est))[5:6] # true estimate = -1.9812346 -0.9895642
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[2:3, 2:3] 

beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
gamma2 = as.matrix(fixef(fit_pbo_est))[5:6] # 1.994051 1.003545
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[2:3, 2:3]

unique_dat = unique(dat_est[,c('subj','w')])
mu1_est = c(); mu2_est = c()
for(i in 1:dim(unique_dat)[1]){
  #if(i %% 10 ==0) print(i)
  mu1 = beta1 + gamma1 * unique_dat$w[i] # mu1 = beta_drg + gamma_drg * unique_dat$w[i]
  mu2 = beta2 + gamma2 * unique_dat$w[i] # mu2 = beta_pbo + gamma_pbo * unique_dat$w[i]
  mu1_est = c(mu1_est, mu1)
  mu2_est = c(mu2_est, mu2)
}
mus = data.frame(mu1_est = mu1_est, mu2_est = mu2_est)

mean(abs(mus$mu1_est - mus$mu2_est)) # 6.156634


load('mc_data_1000_0317.RData')
