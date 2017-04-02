library(Rcpp)
library(foreach)
library(doMC)
library(rbenchmark)
sourceCpp("notebooks/custom_covar.cpp")

x <- matrix(runif(n = 1000 * 1000), nrow = 1000, ncol = 1000)

rbenchmark::benchmark(custom_multi_covar(x), replications = 10)
rbenchmark::benchmark(cov(x), replications = 10)

a <- custom_multi_covar(x)
b <- cov(x)

err.std.cov.against.stable.cov <- max(abs((b - a)/a))

print(err.std.cov.against.stable.cov)

##########

x <- matrix(runif(n = 1000 * 1000), nrow = 1000, ncol = 1000)

distributed.covar <- function(x, n.splits = 16, n.cores = 2) {
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
    s <- custom_multi_covar(x[i:k,])
    rbind(mn, s)
  }
  
  u <- do.call(cbind, res)
  combine_covs(u, ns)
}

rbenchmark::benchmark(distributed.covar(x, n.splits = 4, n.cores = 4), replications = 3)
rbenchmark::benchmark(cov(x), replications = 3)

a <- distributed.covar(x, n.splits = 4, n.cores = 4)
a <- a[2:NROW(a), ]
b <- cov(x)

err.std.cov.against.stable.cov <- max(abs((b - a)/a))

print(err.std.cov.against.stable.cov)

