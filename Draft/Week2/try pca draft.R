### try pca

# continuous variables

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
str(merged)

# continuous variables
con_names = c()
for(i in 1:dim(merged)[2]){
  if(class(merged[, i]) == "integer" | class(merged[, i]) == "numeric"){
    con_names = c(con_names,i)
  }
}
con_data = merged[,con_names]
dim(con_data)
# omit NAs
con_na_data = con_data

for(i in 1:dim(con_data)[2]){
  con_na_data[,i] = ifelse(is.na(con_data[,i]), median(con_data[,i],na.rm = TRUE), con_data[,i])
}

sapply(con_data, function(x) sum(is.na(x)))
sapply(con_na_data, function(x) sum(is.na(x)))


### try pca
dim(con_na_data)
dim(merged)
prin_comp = prcomp(con_na_data,scale. = TRUE)

# cluster the pcs
plot(prin_comp)
plot(prin_comp,type = 'l')
# maybe we can choose the first 4 pc?
dim(prin_comp$x)
pca_data = prin_comp$x[,1:4]
pca_data = as.data.frame(pca_data)
pca_data$ProjectSpecificId = merged$ProjectSpecificId
pca_data$trt = merged$trt
dim(pca_data)
head(pca_data)
km_pca = kmeans(pca_data[,1:4], 2, nstart = 25)
vi(cbind(pca_data$trt, km_pca$cluster))

cluster_summary_table(km_pca,pca_data)
cluster_score_plot2(km_pca,pca_data)

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

cluster_score_plot2 = function(km, data){
  k = length(unique(km$cluster))
  for(clus_num in 1:k){
    id = data[which(km$cluster == clus_num),'ProjectSpecificId']
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
}
