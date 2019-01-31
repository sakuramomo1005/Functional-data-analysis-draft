## 2019-1-31
## try polynomial 
vi = function(x){
  # input: x = nx2 matrix; 1st column first cluster id (1 to k1);
  #            2nd column is the second clustering id (1 to k2)
  n <- dim(x)[1]
  klabels1 <- unique(x[,1])
  klabels2 <- unique(x[,2])
  k1 <- length(klabels1)
  k2 <- length(klabels2)
  # compute the k1 x k2 confusion matrix
  C <- matrix(0,k1,k2)
  for (j1 in klabels1){
    for (j2 in klabels2){
      C[j1,j2] <- sum((x[,1]==j1 & x[,2]==j2)*1)
    }
  }
  p1 <- apply(C,1,sum)/n # compute marginal probabilities
  p2 <- apply(C,2,sum)/n
  H1 <- -sum(na.omit(p1*log(p1)))
  H2 <- -sum(na.omit(p2*log(p2)))
  C <- C/n  
  ICC <- 0
  for (j1 in klabels1){
    for (j2 in klabels2){
      inc <- 0
      if (C[j1,j2] != 0){inc <-C[j1,j2]*log(C[j1,j2]/(p1[j1]*p2[j2]))}
      ICC <- ICC + inc
    }
  }
  vi <- H1 + H2 - 2*ICC
  return(vi)
}
###### Read in data

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week1/first code')

# Read in demographic data
demdat <- read.csv("primary-n151.csv", header=T)
# read in the final hamd outcome data
newdat <- read.csv("longFormat_score17.csv", header=T)

# > dim(demdat)
# [1] 151  67
# > dim(newdat)
# [1] 2009    9

# delete the repeatation of subjects
newdata = newdat[,c(1,2,3,6,7,8,9)]
newdata = unique(newdata)
rownames(newdata) = NULL

dim(newdata) # 287 7 

# merge data
merged = merge(newdata,demdat,by = c('ProjectSpecificId',"site",
                                     "Stage1TX","Sex","age_evaluation"), 
               how = 'inner')
merged$trt = (merged$Stage1TX=="SER/CIT")*1+1
newdat$trt = (newdat$Stage1TX =="SER/CIT")*1+1

dim(merged) # 150 * 7



## fit polynomial models and save the covariates

id = unique(newdat$ProjectSpecificId)
## y ~ x3 + x2 + x1 +1
v = c()
i = 1
data_temp = newdat[newdat$ProjectSpecificId == id[i],]
rownames(data_temp) = NULL
y = data_temp$score17
x = data_temp$week
data = data.frame(x = x, y = y)
ply = lm(y ~ poly(x, 3), data = data)
v = rbind(v,as.vector(ply$coefficients))

new_x = data.frame(x = seq(-10,10,0.01))
#predict(ply,new_x)
plot(new_x[,1], predict(ply,new_x), col = data_temp$trt[1], type = 'l',
     ylim = c(-100,100))
#points(x,y,cex = 0.3)

for(i in 2:length(id)){
  data_temp = newdat[newdat$ProjectSpecificId == id[i],]
  rownames(data_temp) = NULL
  y = data_temp$score17
  x = data_temp$week
  data = data.frame(x = x, y = y)
  ply = lm(y ~ poly(x, 3), data = data)
  v = rbind(v,as.vector(ply$coefficients))
  #print(dim(v))
  #rint(as.vector(ply$coefficients))
  lines(new_x[,1], predict(ply,new_x), col = data_temp$trt[1], type = 'l')
}


## remove nas in the covariates

v = data.frame(v)
v$id = id
v[which(is.na(v[,1])),1] = 0
v[which(is.na(v[,2])),2] = 0
v[which(is.na(v[,3])),3] = 0
v[which(is.na(v[,4])),4] = 0


## kmeans to cluster the covariates
km_plm = kmeans(v[,1:4], 2, nstart = 20)
plotcluster(v[,1:4], km_plm$cluster)

# however, there are very few people in one group

## plot the fitted polynomial lines
y = apply(v[which(km_plm$cluster==1),1:4],2,mean)
new_x = data.frame(x = seq(0,20,0.01))
y = y[1] + new_x* y[2] + new_x^2 * y[3] + new_x^3* y[4] #+ new_x^4 * y[5]+ 
  #new_x^5 * y[6]
plot(new_x[,1],y[,1],col = 1,type ='l', main= 'plot for cluster 1')

y = apply(v[which(km_plm$cluster==2),1:4],2,mean)
new_x = data.frame(x = seq(0,20,1))
y = y[1] + new_x* y[2] + new_x^2 * y[3] + new_x^3* y[4] #+ new_x^4 * y[5]+ 
  #new_x^5 * y[6]
plot(new_x[,1],y[,1],col = 1,type ='l', col = 2, main= 'plot for cluster 2')




## function for the table
cluster_summary_table = function(km, data){
  n_clus = length(unique(km$cluster))
  drg = c(); pbo = c(); tot = c()
  for(i in 1:n_clus){
    drg = c(drg, sum(data[which(km$cluster == i),]$trt == 2,na.rm = TRUE))
    pbo = c(pbo, sum(data[which(km$cluster == i),]$trt == 1,na.rm = TRUE))
    tot = c(tot, dim(data[which(km$cluster == i),])[1])
  }
  drg_pbo = rbind(drg,pbo)
  drg_pbo = rbind(drg_pbo, apply(drg_pbo,2,sum))
  drg_pbo = cbind(drg_pbo, apply(drg_pbo,1,sum))
  drg_pbo = data.frame(drg_pbo)
  rownames(drg_pbo) = c('Drug', 'Placebo', 'Total')
  colnames(drg_pbo) = c(paste('cluster',1:n_clus),'Total')
  return(drg_pbo)
}

cluster_summary_table(km_plm,merged)