library(Rcpp)
library(foreach)
library(doMC)
library(rbenchmark)
sourceCpp("notebooks/custom_mean_var.cpp")

distributed.var <- function(x, n.splits = 16, n.cores = 2) {
  n <- length(x)
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

    custom_mean_var(x[i:k])
  }

  mns <- unlist(lapply(res, function(x) x[1]))
  vars <- unlist(lapply(res, function(x) x[2]))
  return(combine_stats_var(mns, vars, ns))
}

x <- c(runif(10^8, 0, 1))
y <- 0.1 * x + c(runif(10^8, 0, 1))

rbenchmark::benchmark(custom_mean_var(x), replications = 10)
rbenchmark::benchmark(distributed.var(x, n.splits = 32, n.cores = 4), replications = 10)
rbenchmark::benchmark(c(mean(x), sd(x)^2), replications = 10)

a <- custom_mean_var(x)
b <- distributed.var(x, n.splits = 2, n.cores = 2)
c <- c(mean(x), var(x))

err.R.internal.var.wrt.custom.mean.var <- (c - a)/a
err.R.internal.var.wrt.distributed.mean.var <- (c - b)/b
err.distributed.wrt.custom.mean.var <- (b - a)/a

print(err.R.internal.var.wrt.custom.mean.var)
print(err.R.internal.var.wrt.distributed.mean.var)
print(err.distributed.wrt.custom.mean.var)