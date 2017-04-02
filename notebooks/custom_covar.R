library(Rcpp)
library(foreach)
library(doMC)
library(dplyr)

# loading the cpp functions used for numerically stable cov. calculation and 
# also "reduce" operation for the parallel computation of cov.
sourceCpp("notebooks/custom_covar.cpp")

# parallel cov. function
# n.splits (should be of 2^n form) shows the number of splits of data that we want to combine
# n.cores is the number of cores that will be used for the parallel computation
# cov.fun is the base cov. function used on each split. custom_multi_covar is 
#         the 2-pass stable algorithm implemented in custom_covar.cpp
parallel.cov <- function(x, n.splits = 2, n.cores = 2, cov.fun = "custom_multi_covar") {
  x <- x %>% as.matrix()
  n <- NROW(x)
  doMC::registerDoMC(cores = n.cores)
  ns <- rep(round(n/n.splits), n.splits - 1)
  ns <- c(ns, n - sum(ns))
  
  res <- foreach(j = 1:length(ns)) %dopar% {
    if (j == 1) {
      i <- 1
    } else {
      i <- sum(ns[1:(j-1)])+1
    }
    k <- sum(ns[1:j])
    
    mn <- apply(x[i:k,], 2, mean)
    s <- do.call(cov.fun, list(x[i:k,]))
    rbind(mn, s)
  }
  
  u <- do.call(cbind, res)
  mn.cr <- combine_covs(u, ns)
  res <- mn.cr[2:NROW(mn.cr), ]
  rownames(res) <- colnames(x)
  colnames(res) <- colnames(x)
  return(res)
}

# parallel cor. function
# arguments are set similar to parallel.cov
parallel.cor <- function(x, n.splits = 2, n.cores = 2, cov.fun = "custom_multi_covar") {
  cov.mat <- parallel.cov(x, n.splits = n.splits, n.cores = n.cores, cov.fun = cov.fun)
  nrm <- diag(diag(cov.mat)^-0.5)
  cov.mat <- cov.mat %>% as.matrix()
  nrm <- nrm %>% as.matrix()
  res <- nrm %*% cov.mat %*% nrm
  rownames(res) <- colnames(x)
  colnames(res) <- colnames(x)
  return(res)
}
