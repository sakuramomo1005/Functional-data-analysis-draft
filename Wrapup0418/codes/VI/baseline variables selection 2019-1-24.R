### Date 2019-1-24
# Choose baseline variables
# cluster it into k clusters (k = 4 here) and computer the VI with treatment
# labels (drug vs placebo)

# How many baseline variables? 45 continuous
# I tried all combination of one, two, three, and four of them. 


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

dim(merged) # 150 * 7

# how many NAs?
sapply(merged, function(x) sum(is.na(x)))

# this person "MG0060" is not in the 'newdat'
setdiff(demdat$ProjectSpecificId,merged$ProjectSpecificId) # "MG0060"
sum(newdata$ProjectSpecificId == "MG0060")

# choose continuous variables
p = dim(merged)[2] 
merged2 = merged
covars = c()
for(i in 2:p){
  if(class(merged[,i])=="integer" | class(merged[,i]) == "numeric"){
    if(length(unique(merged[,i]))>4){
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
  print(names(merged2)[n_covs])
  names1 = c(names1,names(merged2)[n_covs])
  temp = data.frame(trt = merged2$trt, merged2[,n_covs])
  temp = na.omit(temp)
  km_temp= kmeans(temp[,2:dim(temp)[2]], 4)
  vi1 = c(vi1, vi(cbind(temp$trt, km_temp$cluster)))
}
one_parameter = data.frame(names = names1, vi = vi1)


# two parameters
names2 = c();
vi2 = c()
for(i in 2:(length(covars)-1)){
  for(j in (i+1):length(covars)){
  n_covs = c(i,j)
  print(names(merged2)[n_covs])
  names2 = rbind(names2,names(merged2)[n_covs])
  temp = data.frame(trt = merged2$trt, merged2[,n_covs])
  temp = na.omit(temp)
  km_temp= kmeans(temp[,2:dim(temp)[2]], 4)
  vi2 = c(vi2, vi(cbind(temp$trt, km_temp$cluster)))}
}
two_parameter = data.frame(names1 = names2[,1], names2 = names2[,2], vi = vi2)


# three parameters
names3 = c();
vi3 = c()
for(i in 2:(length(covars)-2)){
  for(j in (i+1):(length(covars)-1)){
    for(l in (j+1):length(covars)){
    n_covs = c(i,j,l)
    #print(names(merged2)[n_covs])
    names3 = rbind(names3,names(merged2)[n_covs])
    temp = data.frame(trt = merged2$trt, merged2[,n_covs])
    temp = na.omit(temp)
    km_temp= kmeans(temp[,2:dim(temp)[2]], 4)
    vi3 = c(vi3, vi(cbind(temp$trt, km_temp$cluster)))}
}}
three_parameter = data.frame(names1 = names3[,1], names2 = names3[,2], 
                             names3 = names3[,3],vi = vi3)


## this may takes more then 10 mins
names4 = c();
vi4 = c()
for(i in 2:(length(covars)-3)){
  for(j in (i+1):(length(covars)-2)){
    for(l in (j+1):(length(covars)-1)){
      for(o in (l+1):length(covars)){
      n_covs = c(i,j,l,o)
      #print(names(merged2)[n_covs])
      names4 = rbind(names4,names(merged2)[n_covs])
      temp = data.frame(trt = merged2$trt, merged2[,n_covs])
      temp = na.omit(temp)
      km_temp= kmeans(temp[,2:dim(temp)[2]], 4)
      vi4 = c(vi4, vi(cbind(temp$trt, km_temp$cluster)))}
  }}}

four_parameter = data.frame(names1 = names4[,1], names2 = names4[,2], 
                             names3 = names4[,3], names4 = names4[,4],vi = vi4)

save(one_parameter, file ='one_parameter.RData')
save(two_parameter, file ='two_parameter.RData')
save(three_parameter, file ='three_parameter.RData')
save(four_parameter, file ='four_parameter.RData')

## min value and histgram
min(vi1);min(vi2);min(vi3);min(vi4)

hist(vi4,main = 'VI with 4 predictors')
hist(vi3,main = 'VI with 3 predictors')
hist(vi2,main = 'VI with 2 predictors')
hist(vi1,main = 'VI with 1 predictors')

names1[which(vi1 == min(vi1))]
names2[which(vi2 == min(vi2)),]
names3[which(vi3 == min(vi3)),]
names4[which(vi4 == min(vi4)),]

# how about add every variable? not good 
temp =  na.omit(merged2)
km_temp= kmeans(temp[,2:dim(temp)[2]], 4)
vi(cbind(temp$trt, km_temp$cluster))


# two parameters combination
names2 = c();
vi2 = c()
for(i in 2:(length(covars)-1)){
  for(j in (i+1):length(covars)){
    n_covs = c(i,j)
    print(names(merged2)[n_covs])
    names2 = rbind(names2,names(merged2)[n_covs])
    temp = data.frame(trt = merged2$trt, merged2[,n_covs])
    temp = na.omit(temp)
    temp[,4] = temp[,1]*5 + temp[,2]
    km_temp= kmeans(temp[,4], 4)
    vi2 = c(vi2, vi(cbind(temp$trt, km_temp$cluster)))}
}
two_parameter2 = data.frame(names1 = names2[,1], names2 = names2[,2], vi = vi2)
