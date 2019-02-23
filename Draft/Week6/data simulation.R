
## 2019-2-22

# Simulate data from mixed-effect model 
# data determined by hcaf depression

# read in data
library(lme4)
setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/Week3/from dr.tarpey')
source("cvxcluster-0513.R")
dat = read.table("hcaf.dat", header=T)
dim(dat) # 3364 7

# generate X
tt = as.matrix(0:6) # pt = the order of time points
ni = length(tt) # 7
X = cbind(matrix(1, length(tt), 1), tt, tt^2)

# simulate data set 
pbo = dat[dat$trt == 0, ] # placebo group
drg = dat[dat$trt == 1, ] # drug group

set.seed(123)
# get the parameter: beta and D
fit_pbo <- lmer(y ~ t1 + I(t1^2) + (t1+I(t1^2)|subj),              
                data = pbo, REML = FALSE)
fit_drg <- lmer(y ~ t1 + I(t1^2) + (t1+I(t1^2)|subj),              
                data = drg, REML = FALSE)

# random effect
sigma.drg <- attr(VarCorr(fit_drg), "sc") 
sigma.pbo <- attr(VarCorr(fit_pbo), "sc") 

# beta
beta_pbo = as.matrix(fixef(fit_pbo)[1:3])
beta_drg = as.matrix(fixef(fit_drg)[1:3])

# generate gamma randomly
Gamma_drg=matrix(c(0,-.5,-.1),3,1)
Gamma_pbo=matrix(c(0,.25,.1),3,1)

# the random effct
Dpbo = VarCorr(fit_pbo)$subj[1:3,1:3]
Ddrg=VarCorr(fit_drg)$subj
epbo=eigen(Dpbo)
edrg=eigen(Ddrg)

# generate 1000 subjects in each group
ndrg=100
npbo=100

#Simulate data

# generate baseline covariates
mux = 1
sigmax=1
xdrg=rnorm(ndrg, mean=mux, sd=sigmax) # baseline covariates
xpbo=rnorm(npbo, mean=mux, sd=sigmax)

# combine dataset
datsim=NULL
for (i in 1:ndrg){
  bi <- edrg$vectors%*%sqrt(abs(diag(edrg$values)))%*%as.matrix(rnorm(3))
  bi = bi 
  yi <- X%*%(beta_drg+bi+Gamma_drg*xdrg[i]) + sigma.drg*rnorm(ni) 
  dati <- NULL
  dati$subj <- paste("drg",rep(i,ni), sep="")
  dati$y <- yi
  dati$tt<- tt
  dati$x <- rep(xdrg[i],ni)
  dati$trt <- rep("drg",ni)
  dati <- data.frame(dati)
  datsim <- rbind(datsim, dati)
}
for (i in 1:npbo){
  bi <- epbo$vectors%*%sqrt(abs(diag(epbo$values)))%*%as.matrix(rnorm(3)) # abs? since there are some eigen values small than 0
  bi = bi
  yi <- X%*%(beta_pbo+bi+Gamma_pbo*xpbo[i]) + sigma.pbo*rnorm(ni)
  dati <- NULL
  dati$subj <- paste("pbo",rep(i,ni), sep="")
  dati$y <- yi
  dati$tt<- tt
  dati$x <- rep(xpbo[i],ni)
  dati$trt <- rep("pbo",ni)
  dati <- data.frame(dati)
  datsim <- rbind(datsim, dati)
}

head(datsim)

# estimate the simulated dataset
pbosim = datsim[datsim$trt == 'pbo', ]
drgsim = datsim[datsim$trt == 'drg', ]
fitdrg.sim = lmer(y ~ tt + I(tt^2) +x + x * tt +
                    x * I(tt^2) + (tt+I(tt^2)|subj),
                  data = drgsim, REML = FALSE)
fitpbo.sim =  lmer(y ~ tt + I(tt^2) + x + x * tt +
                     x * I(tt^2) + (tt+I(tt^2)|subj),
                   data = pbosim, REML = FALSE)
# the estimated values are similar with the known input values 
fixef(fitdrg.sim)
Gamma_drg
beta_drg

fixef(fitpbo.sim)
Gamma_pbo
beta_pbo


plot(pbosim$tt, pbosim$y, col = 1, cex = 0.3)
points(drgsim$tt+0.1, drgsim$y, col = 2, cex = 0.3)

# calculate the purity
# the function for purity
pw = function(w){
  # estimated beta, gamma, sigma parameter
  beta_drg_est = as.matrix(fixef(fitdrg.sim))[2:3]
  gamma_drg_est = as.matrix(fixef(fitdrg.sim))[5:6]
  D_drg_est = as.matrix(VarCorr(fitdrg.sim)$subj)[2:3, 2:3]
  
  beta_pbo_est = as.matrix(fixef(fitpbo.sim))[2:3]
  gamma_pbo_est = as.matrix(fixef(fitpbo.sim))[5:6]
  D_pbo_est = as.matrix(VarCorr(fitpbo.sim)$subj)[2:3, 2:3]
  
  # integrate on zi 
  PDF_fzw = function(X) {
    
    mu_drg = beta_drg_est + gamma_drg_est * w 
    mu_pbo = beta_pbo_est + gamma_pbo_est * w 
    
    Q_drg = (-1/2)*t(X-mu_drg)%*%solve(D_drg_est)%*%(X-mu_drg)
    Q_pbo = (-1/2)*t(X-mu_pbo)%*%solve(D_pbo_est)%*%(X-mu_pbo)
    
    # pdf of multivariable normal distribution 
    f_drg = (1/(2*pi))*(1/sqrt(det(D_drg_est)))*exp(Q_drg)
    f_pbo = (1/(2*pi))*(1/sqrt(det(D_pbo_est)))*exp(Q_pbo)
    
    # purity function
    if((f_drg + f_pbo)!=0){
      res = (f_pbo - f_drg)^2 / (f_drg + f_pbo)
    }else{
      res = 0
    }
    return(res)
  }
  return( adaptIntegrate(PDF_fzw, 
                         lowerLimit= c(-Inf,-Inf),  
                         upperLimit=c(Inf,Inf))$integral)
}

# calculate purity on each w, get the mean value
datsim_unique = unique(data.frame(subj = datsim$subj, trt = datsim$trt, w = datsim$x))
head(datsim_unique)

puritys = c()
for(i in 1:dim(datsim_unique)[1]){
  w = datsim_unique$w[i]
  print(i)
  puritys = c(puritys, pw(w))
}

mean(puritys) #[1] 1.281726
save(puritys,file = 'puritys190222.RData')

# visualization
datsim_unique = unique(data.frame(subj = datsim$subj, trt = datsim$trt, w = datsim$x))

X_t = c(0:6)
par(mfrow = c(1,1))
plot(X_t, rep(0,7),col = 'white', ylim = c(-50,10),
     xlab = 'Time',ylab = 'Score', main = 'Trajectory')
for(i in 1:100){
  w = datsim_unique$w[i]
  mu_drg = coef(fitdrg.sim)$subj[i,2:3] + coef(fitdrg.sim)$subj[i,5:6] * w
  mu_drg = unlist(mu_drg)
  y_temp = X[,2:3] %*% c(mu_drg)
  lines(X_t, y_temp, col = datsim_unique$trt[i], lty = 2)
}
for(i in 101:200){
  ii = i - 100
  w = datsim_unique$w[i]
  mu_pbo = coef(fitpbo.sim)$subj[ii,2:3] + coef(fitpbo.sim)$subj[ii,5:6] * w 
  mu_pbo = unlist(mu_pbo)
  y_temp = X[,2:3] %*% mu_pbo
  lines(X_t, y_temp, col = datsim_unique$trt[i])
}

legend(0, -30, legend=c("drug", "placebo"),
       col=c('black', "red"), lty=2:1, cex=0.8)


#  draw the mean plot
X_t = c(0:6)
mean_beta_drg = c()
mean_beta_pbo = c()
for(i in 1:100){
  w = datsim_unique$w[i]
  mean_beta_drg = rbind(mean_beta_drg, coef(fitdrg.sim)$subj[i,2:3] + coef(fitdrg.sim)$subj[i,5:6] * w )
}

for(i in 101:200){
  w = datsim_unique$w[i]
  ii = i - 100
  mean_beta_pbo = rbind(mean_beta_pbo, coef(fitpbo.sim)$subj[ii,2:3] + coef(fitpbo.sim)$subj[ii,5:6] * w )
}
y1 = X[,2:3] %*% apply(mean_beta_drg,2,mean)
y2 = X[,2:3] %*% apply(mean_beta_pbo,2,mean)
plot(X_t, y1, type = 'l', col = 'black', lty = 2, ylim = c(-50,10),
     xlab = 'Time',ylab = 'Score', main = 'Mean Trajectory')
lines(X_t, y2, col = 'red')
legend(0, -30, c("drug", "placebo"),
       col=c('black', "red"), lty=2:1, cex=0.8)

# How about trajectories estimated without baseline covariates

# what if we fit the model without covariates x 
fitdrg.sim2 = lmer(y ~ tt + I(tt^2) + (tt+I(tt^2)|subj),
                   data = drgsim, REML = FALSE)
fitpbo.sim2 = lmer(y ~ tt + I(tt^2) + (tt+I(tt^2)|subj),
                   data = pbosim, REML = FALSE)


# The trajectories: 
datsim_unique = unique(data.frame(subj = datsim$subj, trt = datsim$trt, w = datsim$x))
# draw in one plot
X_t = c(0:6)
par(mfrow = c(1,1))
plot(X_t, rep(0,7),col = 'white', ylim = c(-50,10),
     xlab = 'Time',ylab = 'Score', main = 'Trajectory estimated without baseline covariates')
for(i in 1:100){
  if(datsim_unique$trt[i] == 'drg'){
    mu_drg = coef(fitdrg.sim2)$subj[i,2:3] 
    mu_drg = unlist(mu_drg)
    y_temp = X[,2:3] %*% c(mu_drg)
    lines(X_t, y_temp, col = 'blue' ,lty = 2)
  }}
for(i in 1:100){
  mu_pbo = coef(fitpbo.sim2)$subj[i,2:3] 
  mu_pbo = unlist(mu_pbo)
  y_temp = X[,2:3] %*% mu_pbo
  lines(X_t, y_temp, col = 'darkgreen')
}

legend(0, -30, legend=c("drug_no_baseline", "placebo_no_baseline"),
       col=c('blue','darkgreen'), lty=2:1, cex=0.8)


# The mean trajectories of the placebo and drug groups are:

mean_beta_drg_est = apply(coef(fitdrg.sim2)$subj[,2:3],2,mean)
mean_beta_pbo_est = apply(coef(fitpbo.sim2)$subj[,2:3],2,mean)
y_drg_est_without = X[,2:3] %*% mean_beta_drg_est
y_pbo_est_without = X[,2:3] %*% mean_beta_pbo_est

plot(X_t, y_drg_est_without, type = 'l', lty = 2, col = 'blue', ylim = c(-50,10),xlab = 'Time',ylab = 'Score', main = 'Mean Trajectory estimated without baseline covariates')
lines(X_t, y_pbo_est_without, col = 'darkgreen')
legend(0, -30, legend=c("drug", "placebo"),
       col=c('blue', 'darkgreen'), lty=2:1, cex=0.8)
