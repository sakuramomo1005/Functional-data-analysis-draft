################################# Boundary calculation #################################
# The boudary took consideration of the linear transformation variable w
## 1. read in data
set.seed(123)
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/Week3/from dr.tarpey')
dat = read.table("hcaf.dat", header=T)
dat = dat[dat$trt != 2, ]; rownames(dat) = NULL
head(dat)

## 2. load functions
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week13')
source('functions_190504.R')
source('cvxcluster-0513.R')

## 3. calculate the combination for max purity
# scale the covariates
dat$age = scale(dat$age)
dat$BaselineCGI = scale(dat$BaselineCGI)

puritys = c()
for(i in seq(0, pi, 0.01)){
  print(i)
  A = c(sin(i), cos(i))
  temp = purity_function(A, varname = c('age','BaselineCGI'), 
                         times = 't1', trt ='trt', trtlevel = c(0,1), 
                         subj ='subj', outcome = 'y', data = dat)
  puritys = c(puritys, sum(temp$purity))
}
i = seq(0, pi, 0.01)[which(puritys == max(puritys))]
AA = c(sin(i), cos(i))

## 4. add the combination term
dat$W = cbind(dat$age, dat$BaselineCGI) %*% matrix(AA,2,1)

## 5. calculate the orthogonal transformation
t = as.matrix(0:6) 
ni = length(t) 
X = cbind(matrix(1, length(t), 1), t, t^2)
Xtpo = X
tbar = mean(t)
Xtpo[, 2] = X[, 2] - tbar
Xtpo[, 3] = (t - tbar)^2 - (ni^2 - 1) / 12
c0 = sqrt(sum(Xtpo[,1]^2))
c1 = sqrt(sum(Xtpo[,2]^2))
c2 = sqrt(sum(Xtpo[,3]^2))
Xtpo[,1] = Xtpo[,1] / c0
Xtpo[,2] = Xtpo[,2] / c1
Xtpo[,3] = Xtpo[,3] / c2
A = matrix(0,3,3) # A = transformation matrix
A[1, 1] = 1 / c0
A[1, 2] = - tbar / c1
A[2, 2] = 1 / c1
A[1, 3] = (tbar^2 - (ni^2 - 1) / 12) / c2
A[2, 3] = -2*tbar / c2
A[3, 3] = 1 / c2

## 6. generate the new dataset with only unique record. 
dat2 = unique(dat[,c('subj','trt','responder','W')])
dim(dat2) # 358 4
rownames(dat2) = NULL
head(dat2)

## 7. fit model and calculate beta and gamma values
dat_est = dat
dat_pbo_est = dat_est[dat_est$trt == 0, ]
dat_drg_est = dat_est[dat_est$trt == 1, ]
fit_drg_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                     W * I(t1^2) + (t1+I(t1^2)|subj),
                   data = dat_drg_est, REML = FALSE)
fit_pbo_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                     W * I(t1^2) + (t1+I(t1^2)|subj),
                   data = dat_pbo_est, REML = FALSE)

# beta, gamma, and D matrix
beta1 = as.matrix(fixef(fit_drg_est))[2:3] # -4.2872329  0.3693021
gamma1 = as.matrix(fixef(fit_drg_est))[5:6] # -0.12896847 -0.01670312
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[1:3, 1:3] 
beta2 = as.matrix(fixef(fit_pbo_est))[2:3] # -4.2557207  0.5003109
gamma2 = as.matrix(fixef(fit_pbo_est))[5:6] # 0.11023236 -0.03155707
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[1:3, 1:3] 

# coefficients for each subject
bis_drg = as.matrix(coef(fit_drg_est)$subj) 
dim(bis_drg) # 196
n_drg = dim(bis_drg)[1]
bis_pbo = as.matrix(coef(fit_pbo_est)$subj) 
dim(bis_pbo) # 162
n_pbo = dim(bis_pbo)[1]
bisall = as.data.frame(rbind(bis_pbo, bis_drg))
bisall$group = c(rep(1, n_pbo),rep(2, n_drg))
bisall$subj = rownames(bisall)
bisall = merge(bisall, dat2, by = 'subj')
bisall$W = bisall$W.y
head(bisall)

## 8. calculate the x matrix after counting the linear combination and orthogonal transformation
bisall$X0 = bisall$`(Intercept)` + bisall$W.x * bisall$W.y
bisall$X1 = bisall$t1 + bisall$`t1:W` * bisall$W.y
bisall$X2 = bisall$`I(t1^2)` + bisall$`I(t1^2):W` * bisall$W.y
bis_trans = cbind(bisall$X0, bisall$X1, bisall$X2) %*% t(solve(A))
data_trans = data.frame(bis_trans[,2:3])
data_trans$group = bisall$group
data_trans$W = bisall$W.y
data_trans$responder = bisall$responder
colnames(data_trans)[1:2] = c("slope", "concavity")
head(data_trans)

## 9. calculate the expected MVN distributions for drug group and placebo group
# get the centor beta
beta1 = as.matrix(fixef(fit_drg_est))[1:3]
gamma1 = as.matrix(fixef(fit_drg_est))[4:6] 
D1 = as.matrix(VarCorr(fit_drg_est)$subj)[1:3, 1:3] 
beta2 = as.matrix(fixef(fit_pbo_est))[1:3]
gamma2 = as.matrix(fixef(fit_pbo_est))[4:6] 
D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[1:3, 1:3]
mu1s = matrix(rep(beta1, each = dim(bisall)[1]),dim(bisall)[1],3) + 
  matrix(bisall$W, dim(bisall)[1],1) %*% matrix(gamma1,1,3)
mu2s = matrix(rep(beta2, each = dim(bisall)[1]),dim(bisall)[1],3) + 
  matrix(bisall$W, dim(bisall)[1],1) %*% matrix(gamma2,1,3)
mu1s_trans = t(solve(A) %*% t(mu1s))
mu2s_trans = t(solve(A) %*% t(mu2s))
D1 = solve(A) %*% D1 %*% t(solve(A))
D2 = solve(A) %*% D2 %*% t(solve(A))
dim(mu1s_trans)
dim(data_trans)

## 10. run the Convexity-Based Clustering Approach algorithm
ns1 = ns2 = 100
pi1 = pi2 = 0.5
# Monte Carlo simulation to calculate the boundary

XSIM0 = c()
for(i in 1:dim(mu1s_trans)[1]){
  miu1 = mu1s_trans[i,2:3]
  miu2 = mu2s_trans[i,2:3]
  cov1 = D1[2:3, 2:3]; cov2 = D2[2:3, 2:3]
  d = 2
  
  f1 = list(miu = as.numeric(miu1), cov = as.matrix(cov1))
  f2 = list(miu = as.numeric(miu2), cov = as.matrix(cov2))
  
  xsim_patient = data.frame(mvrnorm(ns1, f1$miu, f1$cov))
  xsim_patient$group = rep(1, ns1)
  xsim_control = data.frame(mvrnorm(ns2, f2$miu, f2$cov))
  xsim_control$group = rep(2, ns2)
  
  xsim0 = rbind(xsim_patient, xsim_control)
  xsim0$lambda = NA
  
  for (i in 1:(ns1+ns2)){
    xsim0$lambda[i] = lambda(xsim0[i,1:d], d, f1,f2,pi1,pi2)# calculate the lambda of each point
  }
  
  XSIM0 = rbind(XSIM0, xsim0)
}

head(XSIM0)
quantile(XSIM0$lambda)

XSIM0 = XSIM0[order(XSIM0$group),]
rownames(XSIM0) = NULL
ns1 = sum(XSIM0$group == 1)
ns2 = sum(XSIM0$group == 2)

# calculate the boundary of the clusters (the boundary calculated with the linear transformation W)
nby = which(names(XSIM0) == "group")
k = 4
niter = 100
p1 = clustering(XSIM0, k, d, niter, ns1, ns2)
x = p1$xsim #simulated observations
u = p1$bound #threshold of the clusters
u # 0.3599037 0.5348294 0.7373700
boundary = u
# save the boundary 
save(boundary, file='boundary.RData')

## 11. generate the lambda for the real data 
data_trans$lambda = NA
for(i in 1:dim(data_trans)[1]){
  temp = data_trans[i,]
  
  miu1 = matrix(fixef(fit_drg_est)[1:3],1,3) + 
    matrix(temp$W, 1,1) %*% matrix(fixef(fit_drg_est)[4:6],1,3)
  miu2 = matrix(fixef(fit_pbo_est)[1:3],1,3) + 
    matrix(temp$W, 1,1) %*% matrix( fixef(fit_pbo_est)[4:6],1,3)
  miu1 = t(solve(A) %*% t(miu1)); miu1 = miu1[2:3]
  miu2 = t(solve(A) %*% t(miu2)); miu2 = miu2[2:3]
  m1 = rbind(m1, miu1); m2 = rbind(m2, miu2)
  input = temp[,c('slope', 'concavity')]
  
  data_trans$lambda[i] = lambda(input, d = 2, f1 = list(miu=miu1, cov=D1[2:3,2:3]), 
                                f2 = list(miu=miu2, cov=D2[2:3,2:3]), 
                                pi1 = 0.5, pi2 = 0.5)
  
  if (!is.na(data_trans$lambda[i])) {
    data_trans$cluster[i] = 1
    for (j in 1:(k-1)){
      if (data_trans$lambda[i] >= u[j]) data_trans$cluster[i] = j + 1
      # u is the boundary of the clusters
    }      
  }
}

save(data_trans, file = 'data_trans.RData')
plot(density(data_trans$lambda))
table(data_trans$cluster)