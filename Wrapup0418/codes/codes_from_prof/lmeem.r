#####################################################################3
#  FUNCTION:  lmeem
#  Fits a linear mixed effects model using the EM algorithm
#  (based on the equations from McLachlan and Krishnan's book, page 191-194)
#
# This function allows the inclusion of covariates into the model!
#
#
#
#   INPUT:
#      q:     number of covariates
#    dat:     data matrix with N rows where
#             N = number of rows = m1+m2+...+m_n
#                    where m_i = number of observations for the ith subject.
#            column 1: subject id code
#            column 2: response vectors stacked vertically
#            column 3 to (3+q): covariate(s)
#            Remaining columns: Z_i design matrix for random effects.
#  error:    criteria for convergence - convergence is met once the
#             squared difference betweeen parameter estimates from
#             one iteration to the next is less than error.
#
#
#   OUTPUT:
#      beta = estimated vector of regression fixed effects
#         D = estimated covariance matrix for the random effects
#      sigma= estimated standard deviation for the error
#         bjs= predicted random effects (BLUPS)
#      loglike= log-likelihood
#        nit = number of iterations
#
#
############################################################################


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
    { print(j)
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



