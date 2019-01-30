##### 
# Look at is the distribution of 
# drug & placebo treated individuals within each of the clusters 
# Make tabulates
# Make plots of outcome curves for each cluster


#####
# load data
#####

setwd("/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/Week2")

load('one_parameter.RData')
load('two_parameter.RData')
load('three_parameter.RData')
load('four_parameter.RData')

setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week1/first code')

# Read in demographic data and merge data
demdat = read.csv("primary-n151.csv", header=T)
# read in the final hamd outcome data
newdat = read.csv("longFormat_score17.csv", header=T)
newdata = newdat[,c(1,2,3,6,7,8,9)]
newdata = unique(newdata)
merged = merge(newdata,demdat,
               by = c('ProjectSpecificId',"site",
                      "Stage1TX","Sex","age_evaluation"),how = 'inner')
merged$trt = (merged$Stage1TX=="SER/CIT")*1+1
newdat$trt = (newdat$Stage1TX =="SER/CIT")*1+1

#### Functions

## function to find the members in each cluster with smallest VI
vi_selection = function(n, data){
  # n: number of variables
  # data: input data set
  data_temp = unfactorize(data[data$vi == min(data$vi),])
  cov_temp = as.character(data_temp[1,1:n])
  temp = data.frame(ProjectSpecificId = merged$ProjectSpecificId,
                    trt = merged$trt, merged[,cov_temp])
  colnames(temp) = c('ProjectSpecificId','trt',cov_temp)
  temp = na.omit(temp)
  km_temp = kmeans(temp[,3:dim(temp)[2]], 4, nstart = 25)
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

## function for the plots
cluster_score_plot = function(km, data){
  par(mfrow = c(2,2))
  for(clus_num in 1:4){
    id = data[which(km$cluster == clus_num),1]
    newdat_temp = newdat[newdat$ProjectSpecificId %in% id,
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

#### results
km_1 = vi_selection(1,one_parameter)
km_2 = vi_selection(2,two_parameter)
km_3 = vi_selection(3,three_parameter)
km_4 = vi_selection(4,four_parameter)

## Tables
cluster_summary_table(km_1$km,km_1$data)
cluster_summary_table(km_2$km,km_2$data)
cluster_summary_table(km_3$km,km_3$data)
cluster_summary_table(km_4$km,km_4$data)

## Plots
png('one variable selection.png')
cluster_score_plot(km_1$km,km_1$data)
dev.off()

png('two variables selection.png')
cluster_score_plot(km_2$km,km_2$data)
dev.off()

png('three variables selection.png')
cluster_score_plot(km_3$km,km_3$data)
dev.off()

png('four variables selection.png')
cluster_score_plot(km_4$km,km_4$data)
dev.off()
