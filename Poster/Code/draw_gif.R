# codes to draw the gif
# 2019-05-06


library(MASS)
library(mixtools)

setwd('')
source('sim_data_generation.R')

# set parameters:
# generate beta randomly
beta_drg = as.matrix(c(0,0.1,1),3,1)
beta_pbo = as.matrix(c(0,-0.1,0.9),3,1)
# generate gamma randomly
gamma_drg=matrix(c(0, 4,1),3,1)
gamma_pbo=matrix(c(0,-1,1),3,1)
# generate sigma
sigma_drg = 1
sigma_pbo = 1
# generate bi
set.seed(21)
eign1 = diag(c(1,0.4,0.5))
eign2 = matrix(runif(9,0,1),3,3)
a = eign2 %*% eign1 %*% solve(eign2)
bi_sigma = t(a) %*% (a)
eigen(bi_sigma)
bi_sigma2 = bi_sigma

# generate other parameters
tt = as.matrix(0:6) # pt = the order of time points
ni = length(tt) # 7
X = cbind(matrix(1, length(tt), 1), tt, tt^2)
n = 100 # number of subjects in each group  
p = 2 # number of baseline covariates
ni = 7 # number of time points
theta_true = 60/180 * pi # the true theta
alpha = c(sin(theta_true),cos(theta_true))

# generate simulated dataset
data = true_generation(alpha, p, n, ni, tt, X, beta_drg, gamma_drg, bi_sigma,sigma_drg,
                       beta_pbo, gamma_pbo, bi_sigma2,sigma_pbo)
dat_drg = data$dat_drg
dat_pbo = data$dat_pbo
dat = rbind(dat_drg, dat_pbo)
head(dat)

# fit the model and calculated the beta, gamma and D
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

# transformation
Xstart = mvrnorm(10000, rep(0,p), diag(rep(1,p)))

d1 = eigen(D1)
d1 = d1$vectors %*% diag(sqrt(d1$values))
d2 = eigen(D2)
d2 = d2$vectors %*% diag(sqrt(d2$values))

# calculate the plot xlim, and ylim
maxpointsx = c(); maxpointsy = c(); minpointsx = c(); minpointsy = c()
for(i in seq(-1,1,0.2)){
  wvalue = i
  mu1 = beta1 + gamma1 * wvalue
  mu2 = beta2 + gamma2 * wvalue
  
  points1 = Xstart %*% t(d1) 
  points1[,1] = points1[,1] + mu1[1]; points1[,2] = points1[,2] + mu1[2] 
  
  points2 = Xstart %*% t(d2) 
  points2[,1] = points2[,1] + mu2[1]; points2[,2] = points2[,2] + mu1[2]  
  
  maxpointsx = c(maxpointsx,max(points1[,1], points2[,1]))
  maxpointsy = c(maxpointsy,max(points1[,2], points2[,2]))
  minpointsx = c(minpointsx,min(points1[,1], points2[,1]))
  minpointsy = c(minpointsy,min(points1[,2], points2[,2]))
}
xmax = max(maxpointsx)
ymax = max(maxpointsy)
xmin = min(minpointsx)
ymin = min(minpointsy)

m1 = c(); m2 = c()
for( i in seq(-6,6,0.1)){
  mu1 = beta1 + gamma1 * i
  mu2 = beta2 + gamma2 * i
  m1 = rbind(m1, mu1)
  m2 = rbind(m2, mu2)
}

counts = 0
for(i in seq(0,1,0.1)){
  counts = counts + 1
  png(paste(counts,'.png',sep = ''))
  wvalue = i
  
  mu1 = beta1 + gamma1 * wvalue
  mu2 = beta2 + gamma2 * wvalue
  
  points1 = Xstart %*% t(d1) 
  points1[,1] = points1[,1] + mu1[1]; points1[,2] = points1[,2] + mu1[2] 
  
  points2 = Xstart %*% t(d2) 
  points2[,1] = points2[,1] + mu2[1]; points2[,2] = points2[,2] + mu2[2]  
  
  cov(points1)
  cov(points2)
  par(pty="s")
  plot(points1[,1],points1[,2], cex = 0.1, asp=1, main = paste('w = ',i),
       xlim = c(xmin,xmax), 
       ylim = c(ymin,ymax),col = 'orangered', xlab = 'beta1', ylab = 'beta2')
  
  points(points2[,1],points2[,2], cex = 0.1, col = 'yellow')
  ellipse(mu=colMeans(points1), sigma=cov(points1), alpha = .05, npoints = 250, col="orangered4") 
  ellipse(mu=colMeans(points2), sigma=cov(points2), alpha = .05, npoints = 250, col="yellow4") 
  
  x1 = colMeans(points1)
  x2 = (4*eigen(D1)$values[1]*eigen(D1)$vectors[,1]+colMeans(points1)) 
  x3 = (4*eigen(D1)$values[2]*eigen(D1)$vectors[,2]+colMeans(points1))
  arrows(x1[1],x1[2],x2[1],x2[2], col = 'blue', lwd = 2, length = 0.1)
  arrows(x1[1],x1[2],x3[1],x3[2], col = 'blue', lwd = 2, length = 0.1)
  
  x1 = colMeans(points2)
  x2 = (4*eigen(D2)$values[1]*eigen(D2)$vectors[,1]+colMeans(points2)) 
  x3 = (4*eigen(D2)$values[2]*eigen(D2)$vectors[,2]+colMeans(points2))
  arrows(x1[1],x1[2],x2[1],x2[2], col = 'blue', lwd = 2, length = 0.1)
  arrows(x1[1],x1[2],x3[1],x3[2], col = 'blue', lwd = 2, length = 0.1)
  
  lines(m1[,1],m1[,2], lty = 2, col = 'blue')
  lines(m2[,1],m2[,2], lty = 2, col = 'blue')
  
  legend('topright',legend = c('Drug','Placebo'), col = c('orangered', 'yellow'), pch = c(20,20))
  dev.off()
}
