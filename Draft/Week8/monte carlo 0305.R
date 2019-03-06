#### 3-5
# monte carlo simulation to estimate the integral


# 1. functions
### data generation
true_generation = function(alpha){
  # alpha
  set.seed(123)
  alpha = as.matrix(alpha,p,1)
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
    # bi = epbo$vectors%*%sqrt(abs(diag(epbo$values)))%*%as.matrix(rnorm(3))
    bi = randomeffcet %*%as.matrix(rnorm(3))
    yi = X%*%(beta_pbo+bi+gamma_pbo*w[1]) + sigma_pbo*rnorm(ni)
    pbo_temp$y = yi
    dat_pbo = rbind(dat_pbo, as.data.frame(pbo_temp))
  }
  
  dat_drg = c()
  for(i in 1:n){
    drg_temp = NULL
    drg_temp$subj = rep(paste('drg',i,sep=''),ni)
    drg_temp$trt = rep('drg',ni)
    baseline = as.matrix(rnorm(p),p,1)
    x1 = baseline[1]; x2 = baseline[2]
    w = rep(t(alpha) %*% baseline,ni)
    drg_temp$x1 = rep(x1,ni); drg_temp$x2 = rep(x2,ni)
    drg_temp$w = w
    drg_temp$tt = tt
    # bi = epbo$vectors%*%sqrt(abs(diag(epbo$values)))%*%as.matrix(rnorm(3))
    bi = randomeffcet %*%as.matrix(rnorm(3))
    yi = X%*%(beta_drg+bi+gamma_drg*w[1]) + sigma_pbo*rnorm(ni)
    drg_temp$y = yi
    dat_drg = rbind(dat_drg, as.data.frame(drg_temp))
  }
  return(list(dat_drg = dat_drg, dat_pbo = dat_pbo))
}

### functions to calculate the purity
monta_carlo_pdf = function(Xstart, mu1, D1, mu2, D2){
  
  mu1 = matrix(rep(mu1, 1000),1000,2,byrow = TRUE)
  mu2 = matrix(rep(mu2, 1000),1000,2,byrow = TRUE)
  
  Q1 = diag((-1/2)*(Xstart-mu1)%*%solve(D1)%*%t(Xstart-mu1))
  Q2 = diag((-1/2)*(Xstart-mu2)%*%solve(D2)%*%t(Xstart-mu2))
  
  f1 = (1/(2*pi))*(1/sqrt(det(D1)))*exp(Q1)
  f2 = (1/(2*pi))*(1/sqrt(det(D2)))*exp(Q2)

  f3 = (f1 + f2)
  f3 = ifelse(f3 == 0, 1e-4,f3)
  purity = c((f1 - f2)^2 / f3) * 4
  return(purity)
}

purity_calculation = function(dat, beta1, beta2, gamma1, gamma2, D1, D2, p=2){
  unique_dat = unique(dat[,c('subj','w')])
  purity = c()
  for(i in 1:dim(unique_dat)[1]){
    #if(i %% 10 ==0) print(i)
    mu1 = beta1 + gamma1 * unique_dat$w[i]
    mu2 = beta2 + gamma2 * unique_dat$w[i]
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
  
  fit_pbo_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                       w * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_pbo_est, REML = FALSE)
  fit_pbo_est; beta_pbo; gamma_pbo
  fit_drg_est = lmer(y ~ tt + I(tt^2) + w + w * tt +
                       w * I(tt^2) + (tt+I(tt^2)|subj),
                     data = dat_drg_est, REML = FALSE)
  fit_drg_est; beta_drg; gamma_drg
  
  beta1 = as.matrix(fixef(fit_drg_est))[2:3]
  gamma1 = as.matrix(fixef(fit_drg_est))[5:6]
  D1 = as.matrix(VarCorr(fit_drg_est)$subj)[2:3, 2:3]
  
  beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
  gamma2 = as.matrix(fixef(fit_pbo_est))[5:6]
  D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[2:3, 2:3]
  
  b = purity_calculation(dat_est, beta1, beta2, gamma1, gamma2, D1, D2, p=2)
  
  return((b))
}





# set parameters
tt = as.matrix(0:6) # pt = the order of time points
ni = length(tt) # 7
X = cbind(matrix(1, length(tt), 1), tt, tt^2)

n = 100 # number of subjects in each group  
p = 2 # number of baseline covariates

ni = 7

sigma_drg = 3 
sigma_pbo = 4 

# generate beta randomly
beta_drg = as.matrix(c(0,25,1),3,1)
beta_pbo = as.matrix(c(1,-5,-1),3,1)

# generate gamma randomly
gamma_drg=matrix(c(0,-2,-1),3,1)
gamma_pbo=matrix(c(0,2,1),3,1)

randomeffcet = matrix(c(1,1,2,0.1,0.3,0.2,0,0,0.1),3,3)

Xstart = cbind(runif(1000,-1,1),runif(1000,-1,1))

theta_true = 60/180 * pi
data = true_generation(c(2*sin(theta_true),2*cos(theta_true)))
dat_drg = data$dat_drg
dat_pbo = data$dat_pbo
dat = rbind(dat_drg, dat_pbo)
head(dat)
dim(dat)



x = 60
theta_value = x/180*pi
A = c(sin(theta_value),cos(theta_value))
sum(purity_function(A)) # 6.353625


x = 40
theta_value = x/180*pi
A = c(sin(theta_value),cos(theta_value))
sum(purity_function(A)) # 6.128329

x = 70
theta_value = x/180*pi
A = c(sin(theta_value),cos(theta_value))
sum(purity_function(A)) # 6.095003

f = function(x){
  theta_value = x/180*pi
  A = c(sin(theta_value),cos(theta_value))
  return(sum(purity_function(A)))
}

test_value = c()
for(i in 0:360){
  print(i)
  test_value = c(test_value,f(i))
}
data = data.frame(x = 0:360, y = test_value)

plot(0:360, data$y, type = 'l', cex = 0.2)
points(180, f(180))

save(data,file = 'purity_true_mc_called_data_0305.RData')

# true data 

dat$x1 = dat$age; dat$x2 = dat$BaselineCGI; dat$tt = dat$t1
dat$trt = ifelse(dat$trt == 1,'drg','pbo')
