# Variation of Information metric comparing how well two clusterings
# agree with one another (Meila 2007 JMVA)
# One can change the base for the logarithms
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
