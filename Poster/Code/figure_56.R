# figure 5 and figure 6

for(read in 'readdata'){
  # 1. read data and simulate data
  set.seed(123)
  library(lme4)
  
  setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week13')
  source("cvxcluster-0513.R")
  setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/Week3/from dr.tarpey')
  dat = read.table("hcaf.dat", header=T)
  setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week13')
  dat = dat[dat$trt!=2,]; rownames(dat) = NULL
  dat$X1 = scale(dat$age)
  dat$X2 = scale(dat$BaselineCGI)
  AA = c(0,1)
  dat$W = cbind(dat$X1, dat$X2) %*% matrix(AA,2,1)
  
  t = as.matrix(0:6) # pt = the order of time points
  ni = length(t) # 7
  X = cbind(matrix(1, length(t), 1), t, t^2)
  Xtpo = X
  tbar = mean(t) # 3
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
  
  data = unique(dat[,c('subj','trt','responder','W')])
  dim(data) # 358   4
  rownames(data) = NULL
  
  # fit the LMM
  dat_est = dat
  dat_pbo_est = dat_est[dat_est$trt == 0, ]
  dat_drg_est = dat_est[dat_est$trt == 1, ]
  fit_drg_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                       W * I(t1^2) + (t1+I(t1^2)|subj),
                     data = dat_drg_est, REML = FALSE)
  fit_pbo_est = lmer(y ~ t1 + I(t1^2) + W + W * t1 +
                       W * I(t1^2) + (t1+I(t1^2)|subj),
                     data = dat_pbo_est, REML = FALSE)
  
  beta1 = as.matrix(fixef(fit_drg_est))[2:3]
  gamma1 = as.matrix(fixef(fit_drg_est))[5:6] 
  D1 = as.matrix(VarCorr(fit_drg_est)$subj)[1:3, 1:3] 
  beta2 = as.matrix(fixef(fit_pbo_est))[2:3]
  gamma2 = as.matrix(fixef(fit_pbo_est))[5:6]
  D2 = as.matrix(VarCorr(fit_pbo_est)$subj)[1:3, 1:3]

    bis_drg = as.matrix(coef(fit_drg_est)$subj) 
  dim(bis_drg) # 196
  n_drg = dim(bis_drg)[1]
  bis_pbo = as.matrix(coef(fit_pbo_est)$subj) 
  dim(bis_pbo) # 162
  n_pbo = dim(bis_pbo)[1]
  bisall = as.data.frame(rbind(bis_pbo, bis_drg ))
  bisall$group = c(rep(1, n_pbo),rep(2, n_drg))
  bisall$subj = rownames(bisall)
  bisall = merge(bisall, data, by = 'subj')
  head(bisall)
  bisall$W = mean(bisall$W.y)
  
  # made new X0, X1, X2, with combination of baseline variables, which are set as the mean values
  bisall$X0 = bisall$`(Intercept)` + bisall$W.x * bisall$W
  bisall$X1 = bisall$t1 + bisall$`t1:W` * bisall$W
  bisall$X2 = bisall$`I(t1^2)` + bisall$`I(t1^2):W` * bisall$W
  bis_trans = cbind(bisall$X0, bisall$X1, bisall$X2) %*% t(solve(A))
  data_trans_b = data.frame(bis_trans[,2:3])
  data_trans_b$group = bisall$group
  data_trans_b$W = bisall$W
  data_trans_b$responder = bisall$responder
  colnames(data_trans_b)[1:2] = c("slope", "concavity")
  head(data_trans_b)
  #save(data_trans_b, file = 'data_trans_b.RData')
  
  # get the centor beta
  bisall$W = bisall$W
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
  
  # run the cvxclustr code
  dim(mu1s_trans)
  dim(mu2s_trans)
  dim(data_trans_b)
  ns1 = ns2 = 100
  pi1 = pi2 = 0.5
  XSIM_scaled0 = c()
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
    XSIM_scaled0 = rbind(XSIM_scaled0, xsim0)
  }
  head(XSIM_scaled0)
  quantile(XSIM_scaled0$lambda)
  
  XSIM_scaled0 = XSIM_scaled0[order(XSIM_scaled0$group),]
  rownames(XSIM_scaled0) = NULL
  ns1 = sum(XSIM_scaled0$group == 1)
  ns2 = sum(XSIM_scaled0$group == 2)
  
  # calculate the boundary of the clusters:
  nby = which(names(XSIM_scaled0) == "group")
  k = 4
  niter = 100
  p1 = clustering(XSIM_scaled0, k, d, niter, ns1, ns2)
  # p1$bound
  # [1] 0.3611642 0.5341726 0.7361042
  
  x_scaledw = p1$xsim
  head(x_scaledw)
  table(x_scaledw$cluster)
  dim(x_scaledw)
  # 1     2     3     4 
  # 22397 20953 16546 11704
  
}

####################################### Figure 5 ############################################
# draw the boundary
#par(pty = 's')
png('figure5.png')
#pdf('figure5.pdf')
plot(x_scaledw[,1:2], type = "n", xlab = "Slope", ylab = "Concavity",
     xlim = c(-20,10), ylim = c(-5,15),
     main = "Fluoxetine-Treated Subjects - Clustering")

# draw the boundary
for (j in 4:2){
  Bcurvej = NULL
  epsilon = 1
  nc = 25
  Bj = as.matrix(x_scaledw[x_scaledw$cluster == j, 1:2])
  By = sort(Bj[, 2])
  By = seq(By[1], By[length(By)], length.out = nc)
  for (ic in 1:nc){
    Bslice = Bj[(Bj[, 2] > By[ic] - epsilon) & (Bj[, 2] < By[ic] + epsilon),]
    if(length(Bslice) > 2){
      if(dim(Bslice)[1]!=0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[,1] == min(Bslice[, 1]), ])
      }
    }else{
      if(length(Bslice) > 0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[1] == min(Bslice[1]) ])
      }
    }
    
  }
  lines(Bcurvej[, 1], Bcurvej[, 2], lwd = 1.5, col = 'blue')
}

points(apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'orange')

points(apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'red')

text(-15,10,'C1', cex = 1.1)
text(-12,10,'C2', cex = 1.1)
text(-8,10,'C3', cex = 1.1)
text(-5,10,'C4', cex = 1.1)
dev.off()


####################################### Figure 6 ############################################
# draw the boundary
#par(pty = 's')
# figure 6, control points
#png('figure6_control.png')
pdf('figure6_control.pdf')
plot(asp = 1, x_scaledw[,1:2], type = "n", xlab = "Slope", ylab = "Concavity",
     xlim = c(-20,10), ylim = c(-5,15),
     main = "Fluoxetine-Control Subjects - Clustering")
# draw the boundary
for (j in 4:2){
  Bcurvej = NULL
  epsilon = 1
  nc = 25
  Bj = as.matrix(x_scaledw[x_scaledw$cluster == j, 1:2])
  By = sort(Bj[, 2])
  By = seq(By[1], By[length(By)], length.out = nc)
  for (ic in 1:nc){
    Bslice = Bj[(Bj[, 2] > By[ic] - epsilon) & (Bj[, 2] < By[ic] + epsilon),]
    if(length(Bslice) > 2){
      if(dim(Bslice)[1]!=0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[,1] == min(Bslice[, 1]), ])
      }
    }else{
      if(length(Bslice) > 0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[1] == min(Bslice[1]) ])
      }
    }
    
  }
  lines(Bcurvej[, 1], Bcurvej[, 2], lwd = 1.5, col = 'blue')
}
# draw the cluster center
points(apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'orange')
points(apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'red')
text(-15,10,'C1', cex = 1.1)
text(-12,10,'C2', cex = 1.1)
text(-8,10,'C3', cex = 1.1)
text(-5,10,'C4', cex = 1.1)
load('data_trans.RData')
head(data_trans)

cluster4 = data_trans[data_trans$cluster == 1 & data_trans$group == 1, ]
cluster4_r = data_trans[data_trans$cluster == 1 & data_trans$group == 1 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 1, col = 'red')
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 20, col = 'red')
dev.off()

png('figure6_treat.png')
#pdf('figure6_treat.pdf')
plot(asp = 1, x_scaledw[,1:2], type = "n", xlab = "Slope", ylab = "Concavity",
     xlim = c(-20,10), ylim = c(-5,15),
     main = "Fluoxetine-Treated Subjects - Clustering")
# draw the boundary
for (j in 4:2){
  Bcurvej = NULL
  epsilon = 1
  nc = 25
  Bj = as.matrix(x_scaledw[x_scaledw$cluster == j, 1:2])
  By = sort(Bj[, 2])
  By = seq(By[1], By[length(By)], length.out = nc)
  for (ic in 1:nc){
    Bslice = Bj[(Bj[, 2] > By[ic] - epsilon) & (Bj[, 2] < By[ic] + epsilon),]
    if(length(Bslice) > 2){
      if(dim(Bslice)[1]!=0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[,1] == min(Bslice[, 1]), ])
      }
    }else{
      if(length(Bslice) > 0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[1] == min(Bslice[1]) ])
      }
    }
    
  }
  lines(Bcurvej[, 1], Bcurvej[, 2], lwd = 1.5, col = 'blue')
}
# draw the cluster center
points(apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==2,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'orange')
points(apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[1],
       apply(x_scaledw[x_scaledw$cluster==3,1:2],2,mean)[2], cex = 1.5, pch = 19, col = 'red')
text(-15,10,'C1', cex = 1.1)
text(-12,10,'C2', cex = 1.1)
text(-8,10,'C3', cex = 1.1)
text(-5,10,'C4', cex = 1.1)
dev.off()

pdf('allpoints.pdf')
plot(asp = 1, x_scaledw[,1:2], type = "n", xlab = "Slope", ylab = "Concavity",
     xlim = c(-20,10), ylim = c(-5,15),
     main = "Fluoxetine-Subjects - Clustering")
# draw the boundary
for (j in 4:2){
  Bcurvej = NULL
  epsilon = 1
  nc = 25
  Bj = as.matrix(x_scaledw[x_scaledw$cluster == j, 1:2])
  By = sort(Bj[, 2])
  By = seq(By[1], By[length(By)], length.out = nc)
  for (ic in 1:nc){
    Bslice = Bj[(Bj[, 2] > By[ic] - epsilon) & (Bj[, 2] < By[ic] + epsilon),]
    if(length(Bslice) > 2){
      if(dim(Bslice)[1]!=0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[,1] == min(Bslice[, 1]), ])
      }
    }else{
      if(length(Bslice) > 0){
        Bcurvej = rbind(Bcurvej, Bslice[Bslice[1] == min(Bslice[1]) ])
      }
    }
    
  }
  lines(Bcurvej[, 1], Bcurvej[, 2], lwd = 1.5, col = 'blue')
}
# draw the cluster center
cluster4 = data_trans[data_trans$cluster == 4 & data_trans$group == 2, ]
cluster4_r = data_trans[data_trans$cluster == 4 & data_trans$group == 2 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 1, col = 'purple')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 20, col = 'black')
cluster4 = data_trans[data_trans$cluster == 3 & data_trans$group == 2, ]
cluster4_r = data_trans[data_trans$cluster == 3 & data_trans$group == 2 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 1, col = 'blue')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 20, col = 'darkblue')
cluster4 = data_trans[data_trans$cluster == 2 & data_trans$group == 2, ]
cluster4_r = data_trans[data_trans$cluster == 2 & data_trans$group == 2 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 1, col = 'green')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 20, col = 'darkgreen')
cluster4 = data_trans[data_trans$cluster == 1 & data_trans$group == 2, ]
cluster4_r = data_trans[data_trans$cluster == 1 & data_trans$group == 2 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 1, col = 'orange')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 20, col = 'yellow')


# draw the cluster center
cluster4 = data_trans[data_trans$cluster == 4 & data_trans$group == 1, ]
cluster4_r = data_trans[data_trans$cluster == 4 & data_trans$group == 1 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 2, col = 'purple')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 17, col = 'black')
cluster4 = data_trans[data_trans$cluster == 3 & data_trans$group == 1, ]
cluster4_r = data_trans[data_trans$cluster == 3 & data_trans$group == 1 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 2, col = 'blue')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 17, col = 'darkblue')
cluster4 = data_trans[data_trans$cluster == 2 & data_trans$group == 1, ]
cluster4_r = data_trans[data_trans$cluster == 2 & data_trans$group == 1 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 2, col = 'green')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 17, col = 'darkgreen')
cluster4 = data_trans[data_trans$cluster == 1 & data_trans$group == 1, ]
cluster4_r = data_trans[data_trans$cluster == 1 & data_trans$group == 1 & data_trans$responder == 1, ]
points(cluster4[,1],cluster4[,2], cex = 0.8, pch = 2, col = 'red')
points(cluster4_r[,1],cluster4_r[,2], cex = 0.8, pch = 17, col = 'orange')

legend('topright',
       legend = c('treated-cluster 1','treated-cluster 2','treated-cluster 3','treated-cluster 4',
                  'control-cluster 1','control-cluster 2','control-cluster 3',
                  'control-cluster 4'), cex = 0.8,
       col = c('purple','blue','green','orange','purple','blue','green','red'),
       pch = c(1,1,1,1,2,2,2,2))

text(-15,12,'C1', cex = 1.1)
text(-12,12,'C2', cex = 1.1)
text(-8,12,'C3', cex = 1.1)
text(-5,12,'C4', cex = 1.1)
dev.off()

