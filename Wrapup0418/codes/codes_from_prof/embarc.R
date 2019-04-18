# EMBARC data analysis - Fit a cubic B-spline to the longitudinal 0 to 8 week data
# Outcome = HRSD


lmeem <- function(q,dat, error)
{
  
  p <- dim(dat)[2]-q-2
  n  <- length(unique(dat[,1]))
  subj <- as.matrix(unique(dat[,1]))
  
  
  
  # Fit individual curves to each subject to get initial values
  bhat <- NULL
  sigma2 <- NULL
  X <- NULL
  x <- NULL      # collect covariates
  for (i in 1:n){
    prozi <- dat[dat[,1]==subj[i,1],]
    yi <- as.matrix(prozi[,2])
    Zi <- as.matrix(prozi[,(3+q):dim(dat)[2]])
    Xi <- Zi
    if (q>0){
      for (ic in 1:q){
        Xi <- cbind(Xi, prozi[,(2+ic)]*Zi)
      }
    }
    X <- rbind(X, Xi)
    ni <- dim(yi)[1]
    if (dim(Zi)[1] > (dim(Zi)[2]+0)){
      bhati <- solve(t(Zi)%*%Zi)%*%t(Zi)%*%yi
      yhati <- Zi%*%bhati
      ri <- yi-yhati
      sigma2i <- sum(ri^2)/(ni-dim(Zi)[2])
      bhat <- rbind(bhat,t(bhati))
      sigma2 <- rbind(sigma2, c(i,sigma2i))
    }
    if (q>0){
      x <- rbind(x, prozi[1, 3:(2+q)])
    }
  }
  
  # assign initial values to beta, sigma and D
  beta <- rbind(as.matrix(apply(bhat,2,mean)), matrix(0, q*dim(bhat)[2],1))
  D <- cov(bhat)
  sigma <- sqrt(mean(sigma2))
  
  # Begin EM algorithm
  diff=100
  nit <- 0   # record number of iterations
  bjs <- matrix(0,n,p)
  while (diff>error)
  {
    nit <- nit+1
    betaold <- beta
    Dold <- D
    sigmaold <- sigma
    Estep <- NULL   # Store E-step values from equations 5.65 and 5.66
    betabuild1 <- matrix(0,dim(X)[2], dim(X)[2])
    betabuild2 <- matrix(0,dim(X)[2], 1)
    Dtemp <- matrix(0, p,p)
    sigmatemp <- 0
    loglike <- 0
    for (j in 1:n)
    { 
      prozj <- dat[dat[,1]==subj[j,1],]
      yj <- as.matrix(prozj[,2])
      Zj <- as.matrix(prozj[,(3+q):dim(dat)[2]])
      Xj <- X[dat[,1] == subj[j,1],]
      Rj <- diag(dim(Zj)[1])
      nj <- dim(Zj)[1]
      
      if(is.null(dim(Xj))){
        Xj = matrix(Xj, 1, length(Xj))
      }
      #######  E-Step : Equations 5.65 and 5.66 on page 193 of McLachlan and Krishnan  ##################
      bj <- solve(t(Zj)%*%solve(Rj)%*%Zj+sigma^2*solve(D))%*%t(Zj)%*%solve(Rj)%*%(yj-Xj%*%beta)
      bjs[j,] <- t(bj)
      bbtj <- solve(t(Zj)%*%solve(Rj)%*%Zj/sigma^2+solve(D)) + bj%*%t(bj)
      Estep <- rbind(Estep, cbind(matrix(j, dim(Zj)[2],1), bj, bbtj))
      betabuild1 <- betabuild1 + t(Xj)%*%Rj%*%Xj
      betabuild2 <- betabuild2 + t(Xj)%*%solve(Rj)%*%(yj-Zj%*%bj)
      Dtemp <- Dtemp + bbtj
      ej <- yj - Xj%*%beta - Zj%*%bj
      sigmatemp <- sigmatemp +
        sum(diag(t(Zj)%*%solve(Rj)%*%Zj%*%
                   solve(t(Zj)%*%solve(Rj)%*%Zj/sigma^2+solve(D))))+
        t(ej)%*%solve(Rj)%*%ej
      loglike<- loglike + log(det(Zj%*%D%*%t(Zj)+sigma^2*Rj))+
        t(yj-Xj%*%beta)%*%solve(Zj%*%D%*%t(Zj)+sigma^2*Rj)%*%
        (yj-Xj%*%beta)
    }
    loglike <- -0.5*loglike
    
    ####### M-Step : Equations on page 194  #########################################
    
    beta <- solve(betabuild1)%*%betabuild2
    D <- Dtemp/n
    sigma <- as.numeric(sqrt(sigmatemp/dim(dat)[1]))
    
    
    diff=t(beta-betaold)%*%(beta-betaold)+(sigma-sigmaold)^2
    cat("Iteration =", nit, "\n")
    cat("diff = ", diff, "\n")
    cat("log-likelihood = ", loglike, "\n")
    #   cat("beta=",beta, "sigma=",sigma,"D=",D, "\n")
    
    betaold=beta
    sigmaold=sigma
    Dold=D
  }
  
  results <- list("beta" = beta,
                  "D" = D,
                  "sigma " = sigma,
                  "bjs" = bjs,
                  "loglike" = loglike,
                  "nit" = nit)
  results
}


setwd('/Users/yaolanqiu/Desktop/NYU/rotation/Rotation2/week1/first code')
#setwd("~/Box Sync/desktop/nyu/biostat/rotations/spring2019/Kate")
getwd()

`%ni%` <- Negate(`%in%`) # define a "not in" operator

# files for numerical integration
# w <- read.table("GaussLegendre")
# xn <- subset(w$nodes, w$nodes > 0)
# wn = subset(w$weights, w$nodes > 0)
# to use this to integrate a function f(x) from min to max, compute
# sum((max-min)*f((max-min)*xn+min)*wn)

install.packages('refund')
install.packages('vows')
install.packages('fda')

source("lmeem.r")  # load the EM algorithm function
source("vi.r")
require(refund)
require(plyr)
require(mgcv)
library(vows)
library(fda)

# Read in demographic data
demdat <- read.csv("primary-n151.csv", header=T)

# read in the final hamd outcome data
newdat <- read.csv("longFormat_score17.csv", header=T)
names(newdat)
length(unique(newdat$ProjectSpecificId)) # 287

#newdat4 = merge(newdat3, demdat, by.x = 'subject',by.y ="ProjectSpecificId" ,how = 'inner')

#newdat5 = merge(newdat3, demdat, by.x = 'subject',by.y ="ProjectSpecificId" ,how = 'inner')

# delete rows where hamd is missing
newdat <- na.omit(newdat)
# record subjects that have only baseline hamd
subjnull <- NULL
for (isubj in unique(newdat$ProjectSpecificId)){
  dati <- newdat[newdat$ProjectSpecificId == isubj,]
  if (dim(dati)[1]==1){subjnull <- rbind(subjnull, isubj)}
}

# > subjnull
# [,1]    
# isubj "CU0135"
# isubj "MG0125"
# isubj "MG0187"
# isubj "TX0156"
# isubj "TX0178"
# isubj "UM0013"
# isubj "UM0066"
# isubj "UM0092"

newdat <- newdat[newdat$ProjectSpecificId %ni% subjnull,]
#       = newdat[newdat$ProjectSpecificId %in% subjnull == FALSE,]

# rename some of the variables
datnames <- names(newdat)
datnames[1]<- "subject"
datnames[3]<- "TreatmentCode"
datnames[5] <- "hamd"
names(newdat) <- datnames

n <- length(unique(newdat$subject))
n # 279
week <- sort(unique(newdat$week))

# create.bspline.basis: fda package 

# fit a cubic B-spline basis for hamd outcome trajectories
nbasisb2 <- 5  # number of B-spline basis functions to fit hamd trajectories
bspb2 <- create.bspline.basis(range(week), nbasis=nbasisb2, norder=4)

plot(bspb2)

temp = names(newdat)
newdat2 <- NULL  # add B-spline bases to data set 
subjects <- unique(newdat$subject)
kkk= 0
for (isubj in subjects){
  kkk = kkk+1
  print(kkk)
  dati <- newdat[newdat$subject == isubj,]
  dati <- cbind(dati, eval.basis(dati$week,bspb2))
  names(dati) = temp
  newdat2 <- rbind(newdat2, dati)
}

newdat2 <- data.frame(newdat2)
names(newdat2)[((dim(newdat)[2]+1)):dim(newdat2)[2]] <- c("b1", "b2", "b3", "b4", "b5")
names(newdat2)

# dim(newdat2)
# [1] 1769   14
# > dim(newdat)
# [1] 1769    9
                                                                          
# Fit a linear mixed-effects model using the B-spline curves (with no covariates)
error <- 1/10^8
drgdat <- subset(newdat2, TreatmentCode=="SER/CIT")
pbodat <- subset(newdat2, TreatmentCode=="PLA")
dsubjects <- unique(drgdat$subject) # 143
psubjects <- unique(pbodat$subject) # 144
# length(subjects) 287
n1 <- length(unique(dsubjects))
n2 <- length(unique(psubjects))
c(n1,n2)
EMdrg <- lmeem(0, drgdat[,c(1,5,10,11,12,13,14)], error) 
EMpbo <- lmeem(0, pbodat[,c(1,5,10,11,12,13,14)], error) 
names(EMdrg)
#"beta"    "D"       "sigma "  "bjs"     "loglike" "nit"

betahat.drg <- EMdrg$beta  # vector of fixed effects
betahat.drg
bjs.drg <- EMdrg$bjs
EMdrg$D # covariance matrix of random effects

betahat.pbo <- EMpbo$beta  # vector of fixed effects
betahat.pbo
bjs.pbo <- EMpbo$bjs
EMpbo$D # covariance matrix of random effects

# plots for individuals (data points and fitted curves)
tplot=seq(0,8, length.out=100)
drg.beta <- NULL
pdf("cubicEM2.pdf")
plot(dati$week, dati$hamd, ylim=c(0,30),
     xlab="Week", ylab="HRSD",
     main="HRSD vs Week: Cubic B-Spline Fit")
for (i in 1:dim(EMdrg$bjs)[1]){
  dati <- drgdat[drgdat$subject == dsubjects[i],]
  points(dati$week)
  lines(tplot,
        eval.basis(tplot,bspb2)%*%as.matrix(c(betahat.drg[1]+bjs.drg[i,1],
                                              (betahat.drg[2]+bjs.drg[i,2]),
                                              (betahat.drg[3]+bjs.drg[i,3]), 
                                              (betahat.drg[4]+bjs.drg[i,4]),
                                              (betahat.drg[5]+bjs.drg[i,5])),5,1))
  drg.beta <- rbind(drg.beta, c(betahat.drg[1]+bjs.drg[i,1], 
                                betahat.drg[2]+bjs.drg[i,2],
                                betahat.drg[3]+bjs.drg[i,3], 
                                betahat.drg[4]+bjs.drg[i,4],
                                betahat.drg[5]+bjs.drg[i,5]))
  # Sys.sleep(1)
}
dev.off()
drg.beta <- data.frame(drg.beta)
row.names(drg.beta) <- dsubjects

pbo.beta <- NULL
pdf("cubicEMpbo2.pdf")
plot(dati$week, dati$hamd, ylim=c(0,30),
     xlab="Week", ylab="HRSD",
     main="HRSD vs Week: Cubic B-Spline Fit")
for (i in 1:dim(EMpbo$bjs)[1]){
  dati <- pbodat[pbodat$subject == psubjects[i],]
  points(dati$week, dati$hamd, ylim=c(0,30),
       xlab="Week", ylab="HRSD",
       main="HRSD vs Week: Cubic B-Spline Fit")
  lines(tplot,
        eval.basis(tplot,bspb2)%*%as.matrix(c(betahat.pbo[1]+bjs.pbo[i,1],
                                              (betahat.pbo[2]+bjs.pbo[i,2]),
                                              (betahat.pbo[3]+bjs.pbo[i,3]), 
                                              (betahat.pbo[4]+bjs.pbo[i,4]),
                                              (betahat.pbo[5]+bjs.pbo[i,5])),5,1), col=2)
  pbo.beta <- rbind(pbo.beta, c(betahat.pbo[1]+bjs.pbo[i,1], 
                                betahat.pbo[2]+bjs.pbo[i,2],
                                betahat.pbo[3]+bjs.pbo[i,3], 
                                betahat.pbo[4]+bjs.pbo[i,4],
                                betahat.pbo[5]+bjs.pbo[i,5]))
  #  Sys.sleep(1)
}
dev.off()

pbo.beta <- data.frame(pbo.beta)
row.names(pbo.beta) <- psubjects


i=8
#i=2
plot(drgdat$week, drgdat$hamd, ylim=c(0,30), type="n",
     xlab="Week", ylab="HRSD",
     main="HRSD vs Week: Cubic B-Spline Fit")
dati <- drgdat[drgdat$subject == dsubjects[i],]

eval.basis(dati$week,bspb2)

points(dati$week,dati$hamd, pch=19, cex=2.3)
lines(tplot,
      eval.basis(tplot,bspb2)%*%as.numeric(drg.beta[i,]), lwd=3)

# fit a B-spline to individual i separately
Xi <- eval.basis(dati$week,bspb2)
bhati <- solve(t(Xi)%*%Xi)%*%t(Xi)%*%dati$hamd
lines(dati$week, Xi%*%bhati, col=4, lty=2, lwd=2 )

plot(drgdat$week, drgdat$hamd, ylim=c(0,30), type="n",
     xlab="Week", ylab="HRSD",
     main="HRSD vs Week: Cubic B-Spline Fit")
#drgave <- 0*tplot
for (i in 1:dim(EMdrg$bjs)[1]){
  #  dati <- drgdat[drgdat$subject == dsubjects[i],]
  #  points(dati$week,dati$hamd, pch=19)
  lines(tplot,
        eval.basis(tplot,bspb2)%*%as.numeric(drg.beta[i,]))
}


for (i in 1:dim(EMpbo$bjs)[1]){
  lines(tplot,eval.basis(tplot,bspb2)%*%as.numeric(pbo.beta[i,]), col=2, lty=2)
}

drg.ave <- apply(drg.beta,2,mean)
pbo.ave <- apply(pbo.beta,2,mean)
lines(tplot, eval.basis(tplot,bspb2)%*%drg.ave, lwd=5)
lines(tplot, eval.basis(tplot,bspb2)%*%drg.ave, lwd=4, col=4)
lines(tplot, eval.basis(tplot,bspb2)%*%pbo.ave, lwd=5)
lines(tplot, eval.basis(tplot,bspb2)%*%pbo.ave, lwd=4, col="green")

# center the coefficients for HRSD
drg.beta.center <- sweep(drg.beta, 2, apply(drg.beta,2,mean))
pbo.beta.center <- sweep(pbo.beta, 2, apply(pbo.beta,2,mean))




# newdat3 <- NULL
# # reformat data to short form
# for (i in subjects){
#   dati <- newdat2[newdat2$subject==i,]
#   datii <- NULL
#   # fit a B-spline to individual i separately
#   Xi <- eval.basis(dati$week,bspb2)
#   xtx <- (t(Xi)%*%Xi)
#   if (min(eigen(xtx)$value)  > 1/10^15){
#      bhati <- solve(xtx)%*%t(Xi)%*%dati$hamd
#      datii$dat0 <- dati[1,1:9]
#      datii$b1 <- bhati[1]
#      datii$b2 <- bhati[2]
#      datii$b3 <- bhati[3]
#      datii$b4 <- bhati[4]
#      datii$b5 <- bhati[5]
#      datii <- data.frame(datii)
#      names(datii) <- c(names(dati)[1:9], "b1", "b2", "b3", "b4", "b5")
#      newdat3 <- rbind(newdat3, datii)
#   }
# }

# fit quadratic curves
newdat3 <- NULL
# reformat data to short form
for (i in subjects){
  dati <- newdat2[newdat2$subject==i,]
  datii <- NULL
  # fit a B-spline to individual i separately
  Xi <- cbind(1,dati$week,dati$week^2)
  xtx <- (t(Xi)%*%Xi)
  if (min(eigen(xtx)$value)  > 1/10^15){
    bhati <- solve(xtx)%*%t(Xi)%*%dati$hamd
    datii$dat0 <- dati[1,1:9]
    datii$b0 <- bhati[1]
    datii$b1 <- bhati[2]
    datii$b2 <- bhati[3]
    datii <- data.frame(datii)
    names(datii) <- c(names(dati)[1:9], "b0", "b1", "b2")
    newdat3 <- rbind(newdat3, datii)
  }
}

dim(newdat3)
# variation of information
k <- 2 # cluster data into k clusters
km <- kmeans(newdat3[,11:12],k)  # cluster coefficients (don't include intercept)

k=4
km_age = kmeans(newdat3$age_evaluation,k)

# = 
km_age_2 = kmeans(newdat5$age_evaluation.x,k)

bp = data.frame(trt = newdat5$trt,newdat5$bp.w0_1916)
bp = bp[is.na(bp$newdat5.bp.w0_1916)==0,]
km_bp = kmeans(bp$newdat5.bp.w0_1916,k)


# define new treatment variable to take values 1 & 2
newdat3$trt <- (newdat3$TreatmentCode=="SER/CIT")*1+1

vi(cbind(newdat3$trt, km$cluster))


vi(cbind(newdat5$trt, km_age_2$cluster))


vi(cbind(newdat3$trt, km_age$cluster))

vi(cbind(bp$trt, km_bp$cluster))

