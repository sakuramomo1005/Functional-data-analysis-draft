## 2019-04-19
## try 3d , 3 variables in baseline


true_generation = function(alpha, p, n, ni, tt, X, 
                           beta_drg, gamma_drg, bi_sigma,sigma_drg,
                           beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo){
  # alpha
  set.seed(123)
  alpha = as.matrix(alpha,p,1)
  dat_drg = c()
  for(i in 1:n){
    drg_temp = NULL
    drg_temp$subj = rep(paste('drg',i,sep=''),ni)
    drg_temp$trt = rep('drg',ni)
    baseline = as.matrix(rnorm(p,0,1),p,1)
    x1 = baseline[1]; x2 = baseline[2]; x3 =  baseline[3]
    w = rep(t(alpha) %*% baseline,ni)
    drg_temp$x1 = rep(x1,ni); drg_temp$x2 = rep(x2,ni); drg_temp$x3 = rep(x3,ni) 
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
    x1 = baseline[1]; x2 = baseline[2]; x3 =  baseline[3]
    w = rep(t(alpha) %*% baseline,ni)
    pbo_temp$x1 = rep(x1,ni); pbo_temp$x2 = rep(x2,ni); pbo_temp$x3 = rep(x3,ni) 
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


# xiang si xiang xi xiang si xiang si 

# set parameters:
# generate beta randomly
beta_drg = as.matrix(c(0,-1,-3),3,1)
beta_pbo = as.matrix(c(0,3,1),3,1)

# generate gamma randomly
gamma_drg=matrix(c(0,-2,1),3,1)
gamma_pbo=matrix(c(0,1,-2),3,1)

# bi
set.seed(21)
eign1 = diag(c(1,0.4,0.5))
eign2 = matrix(runif(9,0,1),3,3)
a = eign2 %*% eign1 %*% solve(eign2)
bi_sigma = t(a) %*% (a)
eigen(bi_sigma)

# bi
set.seed(7)
eign1 = diag(c(1,0.5,2))
eign2 = matrix(runif(9,0,2),3,3)
bi_sigma2 = eign2 %*%  t(eign2) /10
eigen(bi_sigma2)

# sigma
sigma_drg = 1
sigma_pbo = 1

tt = as.matrix(0:6) # pt = the order of time points
ni = length(tt) # 7
X = cbind(matrix(1, length(tt), 1), tt, tt^2)

n = 50 # number of subjects in each group  
p = 3 # number of baseline covariates
ni = 7 # number of time points

alpha = c(1,1,1)

data = true_generation(alpha, p, n, ni, tt, X, 
                      beta_drg, gamma_drg, bi_sigma,sigma_drg,
                      beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo)
dat_drg = data$dat_drg
dat_pbo = data$dat_pbo
dat = rbind(dat_drg, dat_pbo)
head(dat)

A = c(1,0,1)
res = purity_function(A, varname = c('x1','x2','x3'), times = '', 
                      trt = '',
                      trtlevel = '',
                      subj = '',
                      outcome = 'y',
                      start = 0, data = dat)
sum(res$purity)
res2 = purity_function(alpha, varname = c('x1','x2','x3'), times = '', 
                      trt = '',
                      trtlevel = '',
                      subj = '',
                      outcome = 'y',
                      start = 0, data = dat)
sum(res2$purity)

A = c(3.125,3.152,3.289)
res = purity_function(A, varname = c('x1','x2','x3'), times = '', 
                      trt = '',
                      trtlevel = '',
                      subj = '',
                      outcome = 'y',
                      start = 0, data = dat)
sum(res$purity)



purity_function_for_optim = function(A){
  res = purity_function(A, varname = c('x1','x2','x3'), times = '', 
                        trt = '',
                        trtlevel = '',
                        subj = '',
                        outcome = 'y',
                        start = 0, data = dat)
  return(sum(res$purity))
}

a = Sys.time()
optim( alpha, purity_function_for_optim)
b = Sys.time()
b-a


f = function(A){
  res = purity_function(A, varname = c('x1','x2','x3'), times = '', 
                        trt = '',
                        trtlevel = '',
                        subj = '',
                        outcome = 'y',
                        start = 0, data = dat)
  return(sum(res$purity))
}

result = evolution(n_population, DNA_size, x_bounder, cross_rate, mutate_rate, n_iterations = 200)


theta = pi/3
alpha = c(sin(theta) * sin(theta), sin(theta) * cos(theta), cos(theta))
res0 = purity_function(alpha, varname = c('x1','x2','x3'), times = '', 
                trt = '',
                trtlevel = '',
                subj = '',
                outcome = 'y',
                start = 0, data = dat)

sum(res0$purity)
theta = 0
alpha = c(sin(theta) * sin(theta), sin(theta) * cos(theta), cos(theta)^2)


f = function(theta){
  A = c(sin(theta) * sin(theta), sin(theta) * cos(theta), cos(theta))
  res = purity_function(A, varname = c('x1','x2','x3'), times = '', 
                        trt = '',
                        trtlevel = '',
                        subj = '',
                        outcome = 'y',
                        start = 0, data = dat)
  return(sum(res$purity))
}

optim(pi/3, f, method = 'Brent', lower = 0, upper = 3.14)





f(pi/3)
f(0)
f(1)
f(2.801253)

res = c()
for(i in seq(0,pi,0.01)){
  print(i)
  res = c(res, f(i))
}
a = Sys.time()
optim(pi/3, f)
b = Sys.time()
b-a

a = Sys.time()
optim(pi/3, f, method = 'Brent', lower = 0, upper = 3.15)
b = Sys.time()

plot(seq(0,pi,0.01), res, type = 'l',
     main = 'Purity with theta',
     xlab = 'theta',ylab = 'purity')

points(pi/3, f(pi/3), col = 'red', pch = 18, cex = 0.8)
points(seq(0,pi,0.01)[which(res == max(res))], 
       max(res), pch = 20, cex = 0.8, col = 'blue')
lines(c(pi/3, pi/3), c(0, f(pi/3)), lty = 3, col = 'red')
lines(c(seq(0,pi,0.01)[which(res == max(res))], 
        seq(0,pi,0.01)[which(res == max(res))]), 
      c(0, max(res)), lty = 2, col = 'blue')

legend('topright',  legend=c("True", "Estimated max"),
                          col=c("red", "blue"), lty=2:3, cex=0.8)



lines(c(1.03, 1.03), c(340, 340))
saveres = res


f = function(A){
  #A = c(sin(theta) * sin(theta), sin(theta) * cos(theta), cos(theta)^2)
  res = purity_function(A, varname = c('x1','x2','x3'), times = '', 
                        trt = '',
                        trtlevel = '',
                        subj = '',
                        outcome = 'y',
                        start = 0, data = dat)
  return(sum(res$purity))
}

res2 = c()
for(ii in seq(0,1,0.2)){
  for(jj in seq(0.01,1,0.2)){
    for(kk in seq(0,1,0.2)){
      print(kk)
      res2 = c(res2, f(c(ii,jj,kk)))
    }
  }
}

iis = rep(seq(0,1,0.2), each = 30)
kks = rep(seq(0,1,0.2), 30)
jjs = rep(rep(seq(0,1,0.2), each = 6),5)
trydata = data.frame(i = iis, j = jjs, k = kks, res = res2)
data = trydata

plot(iis, res2)
plot(jjs, res2)
plot(kks, res2)