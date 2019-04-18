r1 = r2 = r3 = c()
for(times in 1:1000){
  set.seed(times)
  phi = rbind(c(2,-1),c(-1,6))
  mu1 = c(0,1); mu2 = c(2,1); mu3 = c(3,3)
  bc1 = mvrnorm(100, mu = mu1, Sigma = phi) 
  bc2 = mvrnorm(100, mu = mu2, Sigma = phi) 
  bc3 = mvrnorm(100, mu = mu3, Sigma = phi) 
  b = rbind(bc1,bc2,bc3)
  B= data.frame(b)
  B$type = c(rep(1,100),rep(2,100),rep(3,100))
  
  t = seq(0,1,0.1)
  clus = c()
  y1 = y2 = y3 = matrix(100*11,100,11)
  for(i in 1:100){
    y1[i,] =  bc1[i,1] + bc1[i,2] * t + rnorm(11,0,0.5)
    y2[i,] =  bc2[i,1] + bc2[i,2] * t + rnorm(11,0,0.5)
    y3[i,] =  bc3[i,1] + bc3[i,2] * t + rnorm(11,0,0.5)
  }
  
  
  x = cbind(rep(1,11),seq(0,1,0.1))
  udv = svd(x)
  u = udv$u; d = diag(udv$d); v = udv$v
  x0 = x %*% v %*% solve(d)
  
  mu = (mu1 + mu2 + mu3)/3
  B = ((mu1 - mu) %*% t(mu1 - mu) + 
         (mu2 - mu) %*% t(mu2 - mu) +
         (mu3 - mu) %*% t(mu3 - mu))/3
  W = phi
  
  library(expm)
  
  w = sqrtm(W)
  
  wbw = solve(w) %*% B %*% solve(w)
  H = eigen(wbw)$vectors
  D = diag(eigen(wbw)$values)
  tau = solve(w) %*% H
  C = diag(c(3.5,1))
  Ctau = C %*% t(tau)
  
  bhat1 = t(solve(t(x0) %*% x0) %*% t(x0) %*% t(y1))
  bhat2 = t(solve(t(x0) %*% x0) %*% t(x0) %*% t(y2))
  bhat3 = t(solve(t(x0) %*% x0) %*% t(x0) %*% t(y3))
  
  bhat21 = t(Ctau %*% t(bc1))
  bhat22 = t(Ctau %*% t(bc2))
  bhat23 = t(Ctau %*% t(bc3))
  
  b = rbind(bc1,bc2,bc3)
  B= data.frame(b)
  B$type = c(rep(1,100),rep(2,100),rep(3,100))
  B$type = as.factor(B$type)
  
  B3 = data.frame(rbind(bhat21,bhat22,bhat23))
  B3$type = B$type
  
  bh1 = t(d %*% v %*% t(bc1))
  bh2 = t(d %*% v %*% t(bc2) )
  bh3 = t(d %*% v %*% t(bc3) )
  
  B22 = data.frame(rbind(bh1,bh2,bh3))
  b22 = data.frame(rbind(bh1,bh2,bh3))
  B22$type = as.factor(B$type)
  
  res1 = kmeans(B[,1:2], centers = 3, nstart = 25)
  res2 = kmeans(B22[,1:2], centers = 3, nstart = 25)
  res3 = kmeans(B3[,1:2], centers = 3, nstart = 25)
  
  true = rbind(mu1,mu2,mu3)
  res1 = res1$centers
  res2 = t(solve(d %*% v) %*% t(res2$centers))
  res3 = t(solve(Ctau) %*% t(res3$centers))
  
  r1 = c(r1, sum(( res1 - true)^2)/3)
  r2 = c(r2, sum(( res2 - true)^2)/3)
  r3 = c(r3, sum(( res3 - true)^2)/3)
}

r = list(r1,r2,r3)
save(r,file = 'r.RData')


setwd("/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week1")
png('figure4.png')
plot(density(r1),main = 'Squared Error Distributions',ylim = c(0,0.2),
     xlab = 'Squared Error', ylab = 'Density',lty = 1)
lines(density(r2),col = 'red',lty=2)
lines(density(r3),col = 'blue',lty = 3)
legend("topright",legend=c('Original Distribution',
                       'Orthogonal Design',
                       'Canonical Transformation'),
       col=c('black',"red", "blue"), box.lwd=0.2, lty = 1:3, cex=0.8)
dev.off()