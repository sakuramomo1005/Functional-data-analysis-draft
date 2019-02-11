library(MASS)

##Estimate the mean and cov of densities
density_estimate = function(data){
  miu = colMeans(data,na.rm = TRUE)
  M = cov(data,use = "pairwise.complete.obs")
  return(list(miu = miu,cov = M))
}

## calculate the lambda 
lambda = function(x, d, f1, f2, pi1, pi2){
  x = as.numeric(x)
  f1x = exp( -d/2 * log(2*pi) - 0.5 * log(abs(det(f1$cov))) - 0.5* t(x - f1$miu) %*% solve(f1$cov) %*% (x-f1$miu) )
  f2x = exp( -d/2 * log(2*pi) - 0.5 * log(abs(det(f2$cov))) - 0.5* t(x - f2$miu) %*% solve(f2$cov) %*% (x-f2$miu) )
  return(pi2 * f2x / (pi1 * f1x + pi2 * f2x))
}

####Simulate Data from the real data
simulate_normal = function(data, group, ns1, ns2) ##d is number of variables, n is number of data simulated per group, pc is data
{
  d = ncol(data)
  pi1 = ns1 / (ns1 + ns2)
  pi2 = ns2 / (ns1 + ns2)
  groupnames = as.matrix(unique(group))
  f1 = density_estimate(data[group == groupnames[1],])
  f2 = density_estimate(data[group == groupnames[2],])  
  xsim_patient = data.frame(mvrnorm(ns1, f1$miu, f1$cov))
  xsim_patient$group = rep(1, ns1)
  xsim_control = data.frame(mvrnorm(ns2, f2$miu, f2$cov))
  xsim_control$group = rep(2, ns2)
  xsim = rbind(xsim_patient, xsim_control)
  xsim$lambda = NA
  for (i in 1:(ns1+ns2)){
    xsim$lambda[i] = lambda(xsim[i, 1:d], d, f1, f2, pi1, pi2)# calculate the lambda of each point
  }
  return(xsim)
}


clustering=function(xsim0, k, d, niter = 100, ns1, ns2){ #k is number of clusters, d is number of variables, xsim is simulated data
  pi1=ns1/(ns1+ns2)
  pi2=ns2/(ns1+ns2)
  xsim=as.data.frame(xsim0)
  km=kmeans(xsim[,1:d], k, iter.max=50) ##initial k-means clustering
  #calculate the stand points according to the clustering km
  w=NULL
  # w are the center of the clusters
  for (j in 1:k){
    w=rbind(w, 
            pi2*mean((km$cluster==j)[(ns1+1):(ns1+ns2)]*1)/
              (pi1*mean((km$cluster==j)[1:ns1]*1)+pi2*mean((km$cluster==j)[(ns1+1):(ns1+ns2)]*1)))
  }
  
  u = rep(0, k-1)  # store cut-points
  purity = rep(NA,niter) # these are two different measurements of the result
  purity2 = rep(NA,niter)
  #start to do the iteration 
  for (iter in 1:niter){
    w = sort(w)
    for (j in 1:(k-1)){
      u[j] = (w[j] + w[j+1]) / 2
    }
    ##re-cluster the dataset according to the cut points
    xsim$cluster = NA
    indx = (xsim$lambda < u[1])*(1:(ns1 + ns2))
    indx = indx[indx > 0]
    xsim$cluster[indx] = 1
    if (k > 2){
      for (j in 2:(k-1)){
        indx = ((xsim$lambda > u[j-1])&(xsim$lambda < u[j]))*(1:(ns1 + ns2))
        indx = indx[indx>0]
        xsim$cluster[indx] = j
      } 
    }
    indx = (xsim$lambda > u[k-1])*(1:(ns1 + ns2))
    indx = indx[indx > 0]
    xsim$cluster[indx] = k
    
    #calculate the stand points according to the clustering
    w = NULL
    purity[iter] = 0
    purity2[iter] = 0
    for (j in 1:k){
      p1j = mean((xsim$cluster==j)[1:ns1]*1)
      p2j = mean((xsim$cluster==j)[(ns1+1):(ns1+ns2)]*1)
      pj = pi1 * p1j + pi2 * p2j
      w = rbind(w, pi2 * p2j / pj)
      purity[iter] = purity[iter] + (pi1 * p1j - pi2 * p2j)^2 / pj
      purity2[iter] = purity2[iter] + (p1j - p2j)^2
    }
    #cat(iter,purity[iter],purity2[iter],"\n")
  }
  head(xsim)
  return(list(xsim = xsim,purity = purity,purity2 = purity2,bound = u))
}



cvxcluster = function(data0 = data.frame(), miu1 = NULL, cov1 = NULL, miu2 = NULL, cov2 = NULL,
                      by = NULL, k, pc = 0, nsim = 10000, niter = 100, ns1 = 10000, ns2 = 10000, ...) UseMethod("cvxcluster")

cvxcluster.default = function(data0 = data.frame(), miu1 = NULL, cov1 = NULL, miu2 = NULL, cov2 = NULL,
                              by = NULL, k, pc = 0, nsim = 10000, niter = 100, ns1 = 10000, ns2 = 10000, ...){
  #nsim is simulated sample size, k is number of clusters
  if (nsim != 10000) {
    ns1 = nsim
    ns2 = nsim
  }
  pi1 = ns1 / (ns1 + ns2)
  pi2 = ns2 / (ns1 + ns2)
  data = as.data.frame(data0)
  
  if (length(data) > 0){
    nby = which(names(data) == by)
    group = data[, nby]
    data = data[,-nby]
    groupnames = as.matrix(unique(group))
    if (nrow(groupnames) > 2) stop("more than two populations are found")
    
    d = ncol(data)
    n = nrow(data)
    if (pc > d) stop("invalid number of principal components, greater than number of variables")
    if (pc > 0) {
      pc0 = princomp(~., data = data, na.action = na.exclude, cor=TRUE)$scores[,1:pc]
      data = pc0
    }
    xsim0 = simulate_normal(data,group,ns1,ns2)
    
    p1 = clustering(xsim0,k,d,niter,ns1,ns2)
    
    #calculate the proportion of real data in clusters
    u = p1$bound
    lambda0 = rep(NA, n)
    cluster0 = rep(NA, n)
    cc = complete.cases(data)
    f1 = density_estimate(data[group == groupnames[1],])
    f2 = density_estimate(data[group == groupnames[2],])  
    for (i in 1:n){
      if (cc[i]) {
        lambda0[i] = lambda(data[i,], d, f1, f2, pi1, pi2)# calculate the lambda of each point
        cluster0[i] = 1
        for (j in 1:(k-1)){
          if (lambda0[i] >= u[j]) cluster0[i] = j + 1
        }
      }
    }
    result = list(n = n, k = k, ns = ns1 + ns2, niter = niter, pc = pc, data = data0, nby = nby, 
                  xsim = p1$xsim, lambda = lambda0, cluster = cluster0, bound = p1$bound, 
                  density = list(density1 = f1, density2 = f2), purity = p1$purity, purity2 = p1$purity2)
  } else{
    
    d = length(as.numeric(miu1))
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
    n = ns1 + ns2
    nby = which(names(xsim0) == "group")
    p1 = clustering(xsim0, k, d, niter, ns1, ns2)
    result = list(n = ns1 + ns2, k = k, ns = ns1 + ns2, niter = niter, pc = pc, data = xsim0, nby = nby, 
                  xsim = p1$xsim, lambda = xsim0$lambda, cluster = p1$xsim$cluster, bound = p1$bound,
                  density = list(density1=f1,density2=f2), purity = p1$purity, purity2 = p1$purity2)
  }  
  result$call = match.call()
  class(result) = "cvxcluster"
  return(result)
}

print.cvxcluster = function(x,...)
{
  cat("Call:\n")
  print(x$call)
  if (x$pc > 0) { 
    cat("\n", x$pc, " principal components used")
  } else {
    cat("\n Original variables used")
  }
  cat("\n",sum(complete.cases(x$data)),"observations used for clustering")
  cat("\n",x$ns,"simulated observations generated")
  cat("\n",x$niter,"iterations done.")
  cat("\n Number of populations in each cluster:\n")
  k=x$k
  s=table( x$data[,x$nby],(x$cluster))
  #row.names(s)=c("Group 1:", "Group 2:")
  for (i in 1:k){
    colnames(s)[i]=paste("Cluster",i)
  }
  print(s)
}




