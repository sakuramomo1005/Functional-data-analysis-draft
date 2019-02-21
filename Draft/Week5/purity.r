
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/Week3/from dr.tarpey')

library(lme4)
source("cvxcluster-0513.R")

dat = read.table("hcaf.dat", header=T)
dim(dat) # 3364 7 

dat$newcov = dat$BaselineCGI

d0 = dat[dat$trt == 0 ,]
d1 = dat[dat$trt == 1 ,]
fit_d0 = lmer(y ~ t1 + I(t1^2) + newcov + newcov * t1 + 
                newcov * I(t1^2) + (t1+I(t1^2)|subj), 
              data = d0, REML = FALSE)
fit_d1 = lmer(y ~ t1 + I(t1^2) + newcov + newcov * t1 + 
                newcov * I(t1^2) + (t1+I(t1^2)|subj),
              data = d1, REML = FALSE)

beta_d0 = as.matrix(fixef(fit_d0)[1:3])
gamma_d0 = as.matrix(fixef(fit_d0)[4:6])
beta_d1 = as.matrix(fixef(fit_d1)[1:3])
gamma_d1 = as.matrix(fixef(fit_d1)[4:6])

beta_d0

gamma_d0

beta_d1

gamma_d1

beta0 = as.matrix(fixef(fit_d0)[1:3])
gamma0 = as.matrix(fixef(fit_d0)[4:6])

beta1 = as.matrix(fixef(fit_d1)[1:3])
gamma1 = as.matrix(fixef(fit_d1)[4:6])

quadratic0 = function(a,b) {
    X = matrix(c(a,b),nrow=2)
    Q = (-1/2)*t(X-mu0)%*%solve(sigma0)%*%(X-mu0)
}
quadratic1 = function(a,b) {
    X = matrix(c(a,b),nrow=2)
    Q = (-1/2)*t(X-mu1)%*%solve(sigma1)%*%(X-mu1)
}
  
PDF = function(x) {
  f0 = (1/(2*pi))*(1/sqrt(det(sigma0)))*exp(quadratic0(x[1],x[2]))
  f1 = (1/(2*pi))*(1/sqrt(det(sigma1)))*exp(quadratic1(x[1],x[2]))
  if((f0 + f1)!=0){
      res = (f0 - f1)^2 / (f0 + f1)
    }else{res = 0}
  return(res)
}

library(cubature)
# get unique w value for each subject
W = data.frame(subj = dat$subj, newcov = dat$newcov)
W = unique(W)

P = c() # save the purity value in P
for(i in W$newcov){
    w = i
    m0 = beta0 + gamma0 * w; m0 = m0[2:3]
    m1 = beta1 + gamma1 * w; m1 = m1[2:3]
    D0 = as.matrix(VarCorr(fit_d0)$subj)[2:3, 2:3]
    D1 = as.matrix(VarCorr(fit_d1)$subj)[2:3, 2:3]
    
    mu0 =  matrix(m0, nrow=2)
    sigma0 = D0
    mu1 = matrix(m1, nrow=2)
    sigma1 = D1
    
    P = c(P,adaptIntegrate(PDF, 
                           lowerLimit= c(-1,-1),  
                           upperLimit=c(1,1))$integral)
  }
mean(P) 
sum(P) 

# create new covariates
# read in data
dat = read.table("hcaf.dat", header=T)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]

cov01 = rnorm(length(unique(d0$subj)),5,1)
cov02 = rnorm(length(unique(d0$subj)),10,1)
newcov0 = data.frame(subj = unique(d0$subj), newcov1 = cov01, newcov2 = cov02)
d0 = merge(d0,newcov0, by = 'subj')

# create new covariates
cov01 = rnorm(length(unique(d1$subj)),10,1)
cov02 = rnorm(length(unique(d1$subj)),5,1)
newcov1 = data.frame(subj = unique(d1$subj), newcov1 = cov01, newcov2 = cov02)
d1 = merge(d1,newcov1, by = 'subj')

# new covariate, which is the combination of the two new covariates
# let's make it a simple summation first
d0$newcov = d0$newcov1 + d0$newcov2 
d1$newcov = d1$newcov1 + d1$newcov2

dat = rbind(d0, d1)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]
head(dat)

# get unique w value for each subject
W = data.frame(subj = dat$subj, newcov = dat$newcov)
W = unique(W)

P = c() # save the purity value in P
for(i in W$newcov){
    w = i
    m0 = beta0 + gamma0 * w; m0 = m0[2:3]
    m1 = beta1 + gamma1 * w; m1 = m1[2:3]
    D0 = as.matrix(VarCorr(fit_d0)$subj)[2:3, 2:3]
    D1 = as.matrix(VarCorr(fit_d1)$subj)[2:3, 2:3]
    
    mu0 =  matrix(m0, nrow=2)
    sigma0 = D0
    mu1 = matrix(m1, nrow=2)
    sigma1 = D1
    
    P = c(P,adaptIntegrate(PDF, 
                           lowerLimit= c(-1,-1),  
                           upperLimit=c(1,1))$integral)
  }
mean(P) 
sum(P)

# create new covariates
# read in data
dat = read.table("hcaf.dat", header=T)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]

cov01 = rnorm(length(unique(d0$subj)),5,1)
cov02 = rnorm(length(unique(d0$subj)),10,1)
newcov0 = data.frame(subj = unique(d0$subj), newcov1 = cov01, newcov2 = cov02)
d0 = merge(d0,newcov0, by = 'subj')

# create new covariates
cov01 = rnorm(length(unique(d1$subj)),10,1)
cov02 = rnorm(length(unique(d1$subj)),5,1)
newcov1 = data.frame(subj = unique(d1$subj), newcov1 = cov01, newcov2 = cov02)
d1 = merge(d1,newcov1, by = 'subj')

# new covariate, which is the combination of the two new covariates
# let's make it a simple summation first
d0$newcov = d0$newcov1 #+ d0$newcov2 
d1$newcov = d1$newcov1 #+ d1$newcov2

dat = rbind(d0, d1)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]
head(dat)

# get unique w value for each subject
W = data.frame(subj = dat$subj, newcov = dat$newcov)
W = unique(W)

P = c() # save the purity value in P
for(i in W$newcov){
    w = i
    m0 = beta0 + gamma0 * w; m0 = m0[2:3]
    m1 = beta1 + gamma1 * w; m1 = m1[2:3]
    D0 = as.matrix(VarCorr(fit_d0)$subj)[2:3, 2:3]
    D1 = as.matrix(VarCorr(fit_d1)$subj)[2:3, 2:3]
    
    mu0 =  matrix(m0, nrow=2)
    sigma0 = D0
    mu1 = matrix(m1, nrow=2)
    sigma1 = D1
    
    P = c(P,adaptIntegrate(PDF, 
                           lowerLimit= c(-1,-1),  
                           upperLimit=c(1,1))$integral)
  }
mean(P) 
sum(P)

# create new covariates
# read in data
dat = read.table("hcaf.dat", header=T)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]

cov01 = rnorm(length(unique(d0$subj)),5,1)
cov02 = rnorm(length(unique(d0$subj)),10,1)
newcov0 = data.frame(subj = unique(d0$subj), newcov1 = cov01, newcov2 = cov02)
d0 = merge(d0,newcov0, by = 'subj')

# create new covariates
cov01 = rnorm(length(unique(d1$subj)),10,1)
cov02 = rnorm(length(unique(d1$subj)),5,1)
newcov1 = data.frame(subj = unique(d1$subj), newcov1 = cov01, newcov2 = cov02)
d1 = merge(d1,newcov1, by = 'subj')

# new covariate, which is the combination of the two new covariates
# let's make it a simple summation first
d0$newcov = 10 * d0$newcov1 #+ d0$newcov2 
d1$newcov = 10 * d1$newcov1 #+ d1$newcov2

dat = rbind(d0, d1)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]
head(dat)

# get unique w value for each subject
W = data.frame(subj = dat$subj, newcov = dat$newcov)
W = unique(W)

P = c() # save the purity value in P
for(i in W$newcov){
    w = i
    m0 = beta0 + gamma0 * w; m0 = m0[2:3]
    m1 = beta1 + gamma1 * w; m1 = m1[2:3]
    D0 = as.matrix(VarCorr(fit_d0)$subj)[2:3, 2:3]
    D1 = as.matrix(VarCorr(fit_d1)$subj)[2:3, 2:3]
    
    mu0 =  matrix(m0, nrow=2)
    sigma0 = D0
    mu1 = matrix(m1, nrow=2)
    sigma1 = D1
    
    P = c(P,adaptIntegrate(PDF, 
                           lowerLimit= c(-1,-1),  
                           upperLimit=c(1,1))$integral)
  }
mean(P) 
sum(P)

# create new covariates
# read in data
dat = read.table("hcaf.dat", header=T)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]

cov01 = rnorm(length(unique(d0$subj)),5,1)
cov02 = rnorm(length(unique(d0$subj)),10,1)
newcov0 = data.frame(subj = unique(d0$subj), newcov1 = cov01, newcov2 = cov02)
d0 = merge(d0,newcov0, by = 'subj')

# create new covariates
cov01 = rnorm(length(unique(d1$subj)),10,1)
cov02 = rnorm(length(unique(d1$subj)),5,1)
newcov1 = data.frame(subj = unique(d1$subj), newcov1 = cov01, newcov2 = cov02)
d1 = merge(d1,newcov1, by = 'subj')

# new covariate, which is the combination of the two new covariates
# let's make it a simple summation first
d0$newcov = d0$newcov1 #+ d0$newcov2 
d1$newcov = 10 * d1$newcov1 #+ d1$newcov2

dat = rbind(d0, d1)
d0 = dat[dat$trt == 0,]
d1 = dat[dat$trt == 1,]
head(dat)

# get unique w value for each subject
W = data.frame(subj = dat$subj, newcov = dat$newcov)
W = unique(W)

P = c() # save the purity value in P
for(i in W$newcov){
    w = i
    m0 = beta0 + gamma0 * w; m0 = m0[2:3]
    m1 = beta1 + gamma1 * w; m1 = m1[2:3]
    D0 = as.matrix(VarCorr(fit_d0)$subj)[2:3, 2:3]
    D1 = as.matrix(VarCorr(fit_d1)$subj)[2:3, 2:3]
    
    mu0 =  matrix(m0, nrow=2)
    sigma0 = D0
    mu1 = matrix(m1, nrow=2)
    sigma1 = D1
    
    P = c(P,adaptIntegrate(PDF, 
                           lowerLimit= c(-1,-1),  
                           upperLimit=c(1,1))$integral)
  }
mean(P) 
sum(P)
