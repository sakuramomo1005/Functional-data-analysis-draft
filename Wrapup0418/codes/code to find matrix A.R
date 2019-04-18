### 2019-2-12
### find the maxtrix A to combine the covariates

setwd('')
library(lme4)
library(splines)
library(fda)  # Use Ramsay's code to obtain design matrices for various
library(mgcv)
source("cvxcluster-0513.R")

dat <- read.table("hcaf.dat", header=T)
dim(dat) # 3364 7 
length(unique(dat$subj)) # 543
length(unique(dat$t1)) # 7

t <- as.matrix(0:6) # pt = the order of time points
ni <- length(t) # 7
X = cbind(matrix(1, length(t), 1), t, t^2)
Xtpo <- X
tbar = mean(t) # 3
Xtpo[, 2] = X[, 2] - tbar
Xtpo[, 3] = (t - tbar)^2 - (ni^2 - 1) / 12
c0 <- sqrt(sum(Xtpo[,1]^2))
c1 <- sqrt(sum(Xtpo[,2]^2))
c2 <- sqrt(sum(Xtpo[,3]^2))
Xtpo[,1] = Xtpo[,1] / c0
Xtpo[,2] = Xtpo[,2] / c1
Xtpo[,3] = Xtpo[,3] / c2
A <- matrix(0,3,3) # A = transformation matrix
A[1, 1] <- 1 / c0
A[1, 2] <- - tbar / c1
A[2, 2] <- 1 / c1
A[1, 3] <- (tbar^2 - (ni^2 - 1) / 12) / c2
A[2, 3] <- -2*tbar / c2
A[3, 3] <- 1 / c2
p <- dim(X)[2]

placebo <- NULL
prozac <- NULL
dat = subset(dat, trt != 2) # We are not using the imi treatment here


### construct the function: purity = f(A)
### the function is based on the hcaf dataset
### only two covariates are considered in the function, age and baseline
f = function(X){
  x1 = X[1]; x2 = X[2]
  dat$covar = x1 * dat$age + x2 * dat$BaselineCGI
  
  
  t <- as.matrix(0:6) # pt = the order of time points
  ni <- length(t) # 7
  X = cbind(matrix(1, length(t), 1), t, t^2)
  Xtpo <- X
  tbar = mean(t) # 3
  Xtpo[, 2] = X[, 2] - tbar
  Xtpo[, 3] = (t - tbar)^2 - (ni^2 - 1) / 12
  c0 <- sqrt(sum(Xtpo[,1]^2))
  c1 <- sqrt(sum(Xtpo[,2]^2))
  c2 <- sqrt(sum(Xtpo[,3]^2))
  Xtpo[,1] = Xtpo[,1] / c0
  Xtpo[,2] = Xtpo[,2] / c1
  Xtpo[,3] = Xtpo[,3] / c2
  A <- matrix(0,3,3) # A = transformation matrix
  A[1, 1] <- 1 / c0
  A[1, 2] <- - tbar / c1
  A[2, 2] <- 1 / c1
  A[1, 3] <- (tbar^2 - (ni^2 - 1) / 12) / c2
  A[2, 3] <- -2*tbar / c2
  A[3, 3] <- 1 / c2
  p <- dim(X)[2]
  
  placebo <- NULL
  prozac <- NULL
  dat = subset(dat, trt != 2) # We are not using the imi treatment here
  
  
  for (jt in unique(dat$trt)){ # fit lme for each arm
    dati <- dat[dat$trt == jt,]
    if(x1 ==0 & x2 ==0){
      fit1 <- lmer(y ~ t1 + I(t1^2) + (t1+I(t1^2)|subj),
                   data = dati, REML = FALSE)
    }else{
      fit1 <- lmer(y ~ t1 + I(t1^2) + (t1+I(t1^2)|subj) +  (t1+I(t1^2)|covar),
                   data = dati, REML = FALSE)
    }
    
    D <- as.matrix(VarCorr(fit1)$subj) # Covariance matrix for random effects
    D <- D[1:p, 1:p]
    beta <- as.matrix(fixef(fit1))
    sigma <- attr(VarCorr(fit1), "sc")
    bis <- as.matrix(coef(fit1)$subj) %*% t(solve(A))
    responder <- NULL # record subjects that are responders or not
    age <- NULL
    BaselineCGI <- NULL
    for (isubj in unique(dati$subj)){
      datisubj <- dati[dati$subj ==isubj,]
      responder <- rbind(responder, datisubj$responder[1])
      age <- rbind(age, datisubj$age[1])
      BaselineCGI <- rbind(BaselineCGI, datisubj$BaselineCGI[1])
    }
    if (jt == 0){
      placebo$n <- length(unique(dati$subj))
      placebo$dat <- dati
      placebo$D <- D
      placebo$beta <- beta
      placebo$sigma <- sigma
      placebo$bis <- bis
      placebo$responder <- responder
      placebo$age <- age
      placebo$BaselineCGI <- BaselineCGI
      placebo$fit <- fit1
    }
    if (jt == 1){
      prozac$n <- length(unique(dati$subj))
      prozac$dat <- dati
      prozac$D <- D
      prozac$beta <- beta
      prozac$sigma <- sigma
      prozac$bis <- bis
      prozac$responder <- responder
      prozac$age <- age
      prozac$BaselineCGI <- BaselineCGI
      prozac$fit <- fit1
    }
  }
  
  beta <- placebo$beta
  D <- solve(A) %*% placebo$D %*% t(solve(A))
  xbar <- solve(A) %*% beta
  pbobeta = as.numeric(xbar[2:3,])
  mu1 = pbobeta
  pboD = D[2:3,2:3]
  
  beta <- prozac$beta
  D <- solve(A) %*% prozac$D %*% t(solve(A))
  xbar <- solve(A) %*% beta
  prozbeta = as.numeric(xbar[2:3, ])
  mu2 = prozbeta
  prozD = D[2:3, 2:3]
  
  data = as.data.frame( rbind(placebo$bis[, 2:3], prozac$bis[, 2:3]) )
  names(data) = c("slope", "concavity")
  data$group = c(rep(1, placebo$n), rep(2, prozac$n))
  
  k = 4 
  p1 = cvxcluster(miu1 = pbobeta, cov1 = pboD, miu2 = prozbeta, cov2 = prozD, 
                  k = 4, nsim = 1000, niter = 20)
  return(p1$purity[20])
}

eps = 0.001

### the first derivative
f_1 = function(X){
  p = length(X)
  y = c()
  for(i in 1:p){
    X1 = X
    X1[i] = X1[i] + eps
    y = c(y, (f(X1) - f(X))/eps)
  }
  print(f(X1))
  return(y)
}

### the partial derivative
f_1_part = function(X,i){
  X1 = X
  X1[i] = X[i] + eps
  y = (f(X1) - f(X))/eps
  print(f(X1))
  return(y)
}

### the second derivative
f_2 = function(X){
  p = length(X)
  y = matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      X1 = X
      X1[j] = X[j] + eps
      f1 = f_1_part(X1,i)  
      f2 = f_1_part(X,i)
      y[i,j] = (f1 - f2)/eps
    }
  }
  return(y)
}

### Newton method
begin = Sys.time()
x_old = c(1,1)
value = c(); x_value = c();x_value2 = c()
set.seed(123)
for(times in 1:100){
  print('********')
  print(times)
  Fx = f_1(x_old)
  value = c(value,Fx)
  x_value = c(x_value, x_old[1])
  x_value2 = c(x_value2, x_old[2])
  F2 = f_2(x_old)
  x_new = x_old - solve(F2) %*% Fx
  if(sum(abs(x_old - x_new)) < 10e-7){
    print('coveraged!')
    break
  }
  x_old = x_new
}
end = Sys.time()
x_old

# 100 times not coverage. Took about 3 hours. 
