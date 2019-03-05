## using genetic algorithm to find the max value in the simulation

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
purity_calculation = function(dat, beta1, beta2, gamma1, gamma2, D1, D2, p=2){
  purity = c()
  for(i in 1:dim(dat)[1]){
    XX = c(dat$tt[i],dat$tt[i]^2)
    mu1 = beta1 + gamma1 * dat$w[i]
    mu2 = beta2 + gamma2 * dat$w[i]
    Q1 = (-1/2)*t(XX-mu1)%*%solve(D1)%*%(XX-mu1)
    Q2 = (-1/2)*t(XX-mu2)%*%solve(D2)%*%(XX-mu2)
    f1 = (1/(2*pi))*(1/sqrt(det(D1)))*exp(Q1)
    f2 = (1/(2*pi))*(1/sqrt(det(D2)))*exp(Q2)
    if((f1 + f2)!=0){
      res = (f1 - f2)^2 / (f1 + f2)
    }else{
      res = 0
    }
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

theta_purity_function = function(x){
  theta_value = x/180*pi
  A = c(sin(theta_value),cos(theta_value))
  return(sum(purity_function(A)))
}

# 2. simulation 
## please do not run, takes a long time
## the results are saved in the RData file

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

# try theta cos and sin


theta_true = 60/180 * pi
data = true_generation(c(2*sin(theta_true),2*cos(theta_true)))
dat_drg = data$dat_drg
dat_pbo = data$dat_pbo
dat = rbind(dat_drg, dat_pbo)
head(dat)
dim(dat)


theta_purity = c()
for(i in 0:360){
  print(i)
  theta_purity = c(theta_purity,theta_purity_function(i))
}

plot(0:360, theta_purity, type = 'l')
points(60,theta_purity_function(60) ,col = 'red',cex = 2, pch = 18)


## use the genetic algorithm to find the max value
source('genetic algorith 0304.R')

n_population = 50
DNA_size = 8
x_bounder = c(0,360)
cross_rate = 0.9 # the common value in this algorithm
mutate_rate = 0.1 # the common value in this algorithm

f = function(x){
  theta_value = x/180*pi
  A = c(sin(theta_value),cos(theta_value))
  return(sum(purity_function(A)))
}

result = evolution(n_population, DNA_size, x_bounder, cross_rate, mutate_rate, n_iterations = 200)

result$x # 226.4062
result$y # 4.120925

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/Week8')

test_value = c()
for(i in 0:360){
  print(i)
  test_value = c(test_value,f(i))
}
data = data.frame(x = 0:360, y = test_value)
save(data,file = 'purity_true_60_called_data_0304.RData')

# plot the purity 
test_value = data$y
result = list(x = 226.4062, y = 4.120925)
plot(0:360, test_value, type = 'l', cex = 0.2)
points(result$x, result$y, col = 'red', pch = 18)
points(result$x -180, result$y, col = 'red', pch = 18)
points(60, f(60), col = 'blue', pch = 17)
points(60 + 180, f(60), col = 'blue', pch = 18)
legend(0,4, c('Estimated max','True value'), col = c('red','blue'), pch = c(18,17),cex= 0.5)
