# Data coding
#  column 1 - patid: patient id 
#  column 2 - invid: site/investigator id 
#  column 3 - protocol: 2 = HCAF 
#  column 4 - arm: 1 = placebo, 2 = fluoxetine varying dose, 3 = imipramine varying dose
#  column 5 - armconv: 0 = placebo, 1 = fluoxetine varying dose, 2 = imipramine varying dose
#  column 6 - visit: visit number (not week) 
#  column 7 - viscod24: visit code, 0 = baseline period, 1 = treatment period 
#  column 8 - dayrelr2: days relative to randomization visit 
#  column 9 - hamdt17: 17 item HamD 
#  column 10 - hamdt21: 21 item HamD 
#  column 11 - cgisev: clinical global impression severity, 1 = normal (not ill), 2 = borderline ill, 
#                    3 = mildly ill, 4 = moderately ill, 5 = markedly ill, 6 = severely ill, 7 = extremely ill 
#  column 12 - cgimp: clinical global impression improvement, 1 = very much improved, 2 = much improved, 
#                    3 = minimally improved, 4 = not improved, 5 = minimally worse, 6 = much worse, 
#                    7 = very much worse 
#  column 13 - age 
#  column 14 - gender: 0 = female, 1 = male 
#  column 15 - ethnic: 0 = caucasian, 1 = black, 2 = hispanic, 3 = asian, 4 = native american, 5 = other 
#  column 16 - height: cm 
#  column 17 - weight: kg 
#  column 18 - BMI 
#  column 19 - basesev: baseline severity from CGI score 
#  column 20 - durcuril: days of current illness 
#  column 21 - numprev: number of previous episodes of psych illness
#  column 22 - h50drsp1: 0 = nonresponder, 1 = responder

vi <- function(x){
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

unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")){
    df[[i]] = as.character(df[[i]])
  } 
  return(df)
}



## function to find the members in each cluster with smallest VI
vi_selection = function(n, data, predata,k){
  # n: number of variables
  # data: input data set
  data_temp = unfactorize(data[data$vi == min(data$vi),])
  cov_temp = as.character(data_temp[1,1:n])
  temp = data.frame(ProjectSpecificId = predata$PATID,
                    trt = predata$ARM, predata[,cov_temp])
  colnames(temp) = c('ProjectSpecificId','trt',cov_temp)
  temp = na.omit(temp)
  km_temp = kmeans(temp[,3:dim(temp)[2]], k, nstart = 25)
  return(list(km = km_temp, data = temp))
}

## function for the table
cluster_summary_table = function(km, data,predata){
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
cluster_score_plot = function(km, data,predata){
  nk = length(unique(km$cluster))
  par(mfrow = c(2,2))
  for(clus_num in 1:nk){
    id = data[which(km$cluster == clus_num),1]
    newdat_temp = predata[predata$PATID %in% id,
                         c('PATID','DAYRELR2','HAMDT17','ARM')]
    plot_temp = newdat_temp[newdat_temp$PATID == id[1],]
    plot(plot_temp$DAYRELR2,
         plot_temp$HAMDT17,
         type = 'l', xlab = 'week',ylab = 'score',
         ylim =c(0, max(predata$HAMDT17,na.rm =TRUE)),
         col = plot_temp$ARM[1], # black: pbo; red: drg 
         main = paste('cluster',clus_num,'score plot'))
    
    for(i in 2:length(id)){
      plot_temp = newdat_temp[newdat_temp$PATID == id[i],]
      lines(col = plot_temp$ARM[1],
            plot_temp$DAYRELR2,
            plot_temp$HAMDT17)
    }
  }
  par(mfrow = c(1,1))
}


setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/data')

dat <- read.table("hcaf.txt",header=T)
dim(dat) # 4616   22

# Kate - perhaps use column 8 as the time variable but delete the
# pre-randomization times (negative values)

dat = dat[dat$DAYRELR2>=0,]
rownames(dat) = NULL
dim(dat) # 3890   22


##############

#####plot y 

##############
id = dat$PATID
i=1
dat_temp = dat[dat$PATID == id[i],]
x = dat_temp$DAYRELR2
y = dat_temp$HAMDT17
plot(x,y,col = dat_temp$ARM[1],type = 'l',
     ylim = c(min(dat$HAMDT17,na.rm=TRUE),max(dat$HAMDT17,na.rm = TRUE)))

for(i in 2:length(id)){
  dat_temp = dat[dat$PATID == id[i],]
  x = dat_temp$DAYRELR2
  y = dat_temp$HAMDT17
  lines(x,y,col = dat_temp$ARM[1])
}

g1 = dat[dat$ARM == 1,c('PATID','DAYRELR2','HAMDT17')]
g2 = dat[dat$ARM == 2,c('PATID','DAYRELR2','HAMDT17')]
g3 = dat[dat$ARM == 3,c('PATID','DAYRELR2','HAMDT17')]

gg1 = reshape(g1, idvar = "PATID", timevar = "DAYRELR2", direction = "wide")
gg2 = reshape(g2, idvar = "PATID", timevar = "DAYRELR2", direction = "wide")
gg3 = reshape(g3, idvar = "PATID", timevar = "DAYRELR2", direction = "wide")

y1 = apply(gg1,2,mean,na.rm=TRUE)[2:dim(gg1)[2]]
y2 = apply(gg2,2,mean,na.rm=TRUE)[2:dim(gg1)[2]]
y3 = apply(gg3,2,mean,na.rm=TRUE)[2:dim(gg1)[2]]

x1 = as.numeric(unlist(strsplit(names(gg1),"\\."))[seq(3,103,2)])
x2 = as.numeric(unlist(strsplit(names(gg2),"\\."))[seq(3,103,2)])
x3 = as.numeric(unlist(strsplit(names(gg3),"\\."))[seq(3,103,2)])

d1 = data.frame(x1 = x1, y1 = y1)
d1 = d1[order(d1$x1),]
d2 = data.frame(x2 = x2, y2 = y2)
d2 = d2[order(d2$x2),]
d3 = data.frame(x3 = x3, y3 = y3)
d3 = d3[order(d3$x3),]

plot(d1$x1, d1$y1, col = 1, type = 'l', ylim = c(0,40))
lines(d2$x2, d2$y2, col = 2)
lines(d3$x3, d3$y3, col = 3)


##############

##### k means covariances

##############

dim(dat)
dat$PATID2 = as.character(dat$PATID) 

dat2 = aggregate(dat[,3:22],  
                 by = list(dat$PATID2), function(x) c(mean=mean(x,na.rm=TRUE)))
names(dat2)[1] ='PATID'
head(dat2)


### choose k = 3 here since ARMs = 3
covars = names(dat2)[5:21]
data = dat2
names1 = c()
vi1 = c()
for(i in 5:length(covars)){
  n_covs = i
  print(names(data)[n_covs])
  temp = data.frame(trt = dat2$ARM, dat2[,n_covs])
  temp = na.omit(temp)
  if(length(unique(dat2[,n_covs]))>2){
    names1 = c(names1,names(data)[n_covs])
    km_temp= kmeans(temp[,2:dim(temp)[2]], 3)
    vi1 = c(vi1, vi(cbind(temp$trt, km_temp$cluster)))
  }
}
one_parameter = data.frame(names = names1, vi = vi1)


covars = names(dat2)[5:21]
data = dat2
names2 = c();
vi2 = c()
for(i in 5:(length(covars)-1)){
  for(j in (i+1):length(covars)){
    n_covs = c(i,j)
    print(names(data)[n_covs])
    names2 = rbind(names2,names(data)[n_covs])
    temp = data.frame(trt = data$ARM, data[,n_covs])
    temp = na.omit(temp)
    km_temp= kmeans(temp[,2:dim(temp)[2]], 3)
    vi2 = c(vi2, vi(cbind(temp$trt, km_temp$cluster)))}
}
two_parameter = data.frame(names1 = names2[,1], names2 = names2[,2], vi = vi2)


#### results

##############

##### summary table and plots

##############
km_1 = vi_selection(1,one_parameter, dat2,k=3)

cluster_summary_table(km_1$km,km_1$data)
cluster_score_plot(km_1$km,km_1$data)