## 2019-01-28
## how to choose k value?

### since the range vi can change based on different k,
### how to make them comparable?


### standardize vi value
### the formular in the paper: vi/(2(log(k))) 

library(cluster)
library(fpc)

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

# how many NAs?
sapply(merged, function(x) sum(is.na(x)))


#### deal with k

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week2')

## function to find the members in each cluster with smallest VI
vi_selection = function(n, data, k){
  # n: number of variables
  # data: input data set
  data_temp = unfactorize(data[data$vi == min(data$vi),])
  cov_temp = as.character(data_temp[1,1:n])
  temp = data.frame(ProjectSpecificId = merged$ProjectSpecificId,
                    trt = merged$trt, merged[,cov_temp])
  colnames(temp) = c('ProjectSpecificId','trt',cov_temp)
  temp = na.omit(temp)
  km_temp = kmeans(temp[,3:dim(temp)[2]], k, nstart = 25)
  return(list(km = km_temp, data = temp))
}

## function for the table
cluster_summary_table = function(km, data){
  n_clus = length(unique(km$cluster))
  drg = c(); pbo = c(); tot = c()
  for(i in 1:n_clus){
    drg = c(drg, sum(data[which(km$cluster == i),]$trt == 2))
    pbo = c(pbo, sum(data[which(km$cluster == i),]$trt == 1))
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

cluster_score_plot = function(km, data,k){
  par(mfrow = c(2,(round(k/2)+1)))
  for(clus_num in 1:k){
    id = data[which(km$cluster == clus_num),1]
    newdat_temp = newdat[which(newdat$ProjectSpecificId %in% id),
                         c('ProjectSpecificId','week','score17','trt')]
    plot_temp = newdat_temp[newdat_temp$ProjectSpecificId == id[1],]
    plot(plot_temp$week,
         plot_temp$score17,
         type = 'l', xlab = 'week',ylab = 'score',
         ylim =c(0, max(newdat$score17,na.rm =TRUE)),
         col = plot_temp$trt[1], # black: pbo; red: drg 
         main = paste('cluster',clus_num,'score plot'))
    
    for(i in 2:length(id)){
      plot_temp = newdat_temp[newdat_temp$ProjectSpecificId == id[i],]
      lines(col = plot_temp$trt[1],
            plot_temp$week,
            plot_temp$score17)
    }
  }
  par(mfrow = c(1,1))
}



#### Do not adjust VI

min_vi1 = c()
### one parameter
for(k in 2:10){
  # choose continuous variables
  print('******')
  print(paste('k = ',k))
  p = dim(merged)[2] 
  merged2 = merged
  covars = c()
  for(i in 2:p){
    if(class(merged[,i])=="integer" | class(merged[,i]) == "numeric"){
      if(length(unique(merged[,i]))>k){
        covars = c(covars, i)
      }
    }
  }
  merged2 = merged[,c(1,covars)]
  merged2$trt = merged$trt
  
  dim(merged2) # 150  47
  
  ## put into k-means and VI calculation
  # one parameters
  names1 = c();
  vi1 = c()
  for(i in 2:length(covars)){
    n_covs = i
    # print(names(merged2)[n_covs])
    names1 = c(names1,names(merged2)[n_covs])
    temp = data.frame(trt = merged2$trt, merged2[,n_covs])
    temp = na.omit(temp)
    km_temp= kmeans(temp[,2:dim(temp)[2]], k, nstart = 20)
    vi1 = c(vi1, vi(cbind(temp$trt, km_temp$cluster)))
  }
  one_parameter = data.frame(names = names1, vi = vi1)
  print("**** min vi ****")
  print(min(one_parameter$vi))
  min_vi1 = c(min_vi1,min(one_parameter$vi))
  res = vi_selection(1,one_parameter,k)
  print(cluster_summary_table(res$km,res$data))
  print("*** p value ***")
  print(fisher.test(cluster_summary_table(res$km,res$data)[1:2,1:k])$p.value)
}


min_vi2 = c()
### two parameter
for(k in 2:10){
  # choose continuous variables
  print('******')
  print(paste('k = ',k))
  p = dim(merged)[2] 
  merged2 = merged
  covars = c()
  for(i in 2:p){
    if(class(merged[,i])=="integer" | class(merged[,i]) == "numeric"){
      if(length(unique(merged[,i]))>k){
        covars = c(covars, i)
      }
    }
  }
  merged2 = merged[,c(1,covars)]
  merged2$trt = merged$trt
  
  dim(merged2) # 150  47
  
  ## put into k-means and VI calculation
  
  names2 = c();
  vi2 = c()
  for(i in 2:(length(covars)-1)){
    for(j in (i+1):length(covars)){
      n_covs = c(i,j)
      #print(names(merged2)[n_covs])
      names2 = rbind(names2,names(merged2)[n_covs])
      temp = data.frame(trt = merged2$trt, merged2[,n_covs])
      #print("****")
      #print(dim(temp))
      temp = na.omit(temp)
      #print(dim(temp))
      km_temp= kmeans(temp[,2:dim(temp)[2]], k, nstart = 25)
      vi2 = c(vi2, vi(cbind(temp$trt, km_temp$cluster)))
    }
  }
  two_parameter = data.frame(names1 = names2[,1], names2 = names2[,2], vi = vi2)
  
  print("**** min vi ****")
  print(min(two_parameter$vi))
  min_vi2 = c(min_vi2,min(two_parameter$vi))
  res = vi_selection(1,two_parameter,k)
  print(cluster_summary_table(res$km,res$data))
  print("*** p value ***")
  print(fisher.test(cluster_summary_table(res$km,res$data)[1:2,1:k])$p.value)
}


#### Adjusted VI

min_vi1_adj = c()
### one parameter
for(k in 2:10){
  # choose continuous variables
  print('******')
  print(paste('k = ',k))
  p = dim(merged)[2] 
  merged2 = merged
  covars = c()
  for(i in 2:p){
    if(class(merged[,i])=="integer" | class(merged[,i]) == "numeric"){
      if(length(unique(merged[,i]))>k){
        covars = c(covars, i)
      }
    }
  }
  merged2 = merged[,c(1,covars)]
  merged2$trt = merged$trt
  
  dim(merged2) # 150  47
  
  ## put into k-means and VI calculation
  # one parameters
  names1 = c();
  vi1 = c()
  for(i in 2:length(covars)){
    n_covs = i
    # print(names(merged2)[n_covs])
    names1 = c(names1,names(merged2)[n_covs])
    temp = data.frame(trt = merged2$trt, merged2[,n_covs])
    temp = na.omit(temp)
    km_temp= kmeans(temp[,2:dim(temp)[2]], k, nstart = 20)
    vi1 = c(vi1, vi(cbind(temp$trt, km_temp$cluster))/(2*log(max(k,2))))
  }
  one_parameter = data.frame(names = names1, vi = vi1)
  print("**** min vi ****")
  print(min(one_parameter$vi))
  min_vi1_adj = c(min_vi1_adj,min(one_parameter$vi))
  res = vi_selection(1,one_parameter,k)
  print(cluster_summary_table(res$km,res$data))
  print("*** p value ***")
  print(fisher.test(cluster_summary_table(res$km,res$data)[1:2,1:k])$p.value)
}


min_vi2_adj = c()
### two parameter
for(k in 2:10){
  # choose continuous variables
  print('******')
  print(paste('k = ',k))
  p = dim(merged)[2] 
  merged2 = merged
  covars = c()
  for(i in 2:p){
    if(class(merged[,i])=="integer" | class(merged[,i]) == "numeric"){
      if(length(unique(merged[,i]))>k){
        covars = c(covars, i)
      }
    }
  }
  merged2 = merged[,c(1,covars)]
  merged2$trt = merged$trt
  
  dim(merged2) # 150  47
  
  ## put into k-means and VI calculation
  
  names2 = c();
  vi2 = c()
  for(i in 2:(length(covars)-1)){
    for(j in (i+1):length(covars)){
      n_covs = c(i,j)
      #print(names(merged2)[n_covs])
      names2 = rbind(names2,names(merged2)[n_covs])
      temp = data.frame(trt = merged2$trt, merged2[,n_covs])
      #print("****")
      #print(dim(temp))
      temp = na.omit(temp)
      #print(dim(temp))
      km_temp= kmeans(temp[,2:dim(temp)[2]], k, nstart = 25)
      vi2 = c(vi2, vi(cbind(temp$trt, km_temp$cluster))/(2*log(max(k,2))))
      }
  }
  two_parameter = data.frame(names1 = names2[,1], names2 = names2[,2], vi = vi2)
  
  print("**** min vi ****")
  print(min(two_parameter$vi))
  min_vi2_adj = c(min_vi2_adj,min(two_parameter$vi))
  res = vi_selection(1,two_parameter,k)
  print(cluster_summary_table(res$km,res$data))
  print("*** p value ***")
  print(fisher.test(cluster_summary_table(res$km,res$data)[1:2,1:k])$p.value)
}
