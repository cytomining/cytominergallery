#' Parallel covariance function
#'
#' @param x ...
#' @param splits shows the number of splits of data that we want to combine  (should be of 2^n form)
#' @param cores number of cores that will be used for the parallel computation
#' @param cov_fun covariance function used on each split
#'
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom foreach %dopar%
#'
parallel_cov <- function(x, splits = 2, cores = 2, cov_fun = "two_pass_multi_covar") {

  x <- as.matrix(x)

  n <- NROW(x)

  doParallel::registerDoParallel(cores = cores)
  
  ns <- rep(round(n/splits), splits - 1)

  ns <- c(ns, n - sum(ns))

  j <- 0 # to avoid warning: no visible binding for global variable ‘j’
  
  result <- foreach::foreach(j = 1:length(ns)) %dopar% {
    if (j == 1) {
      i <- 1

    } else {
      i <- sum(ns[1:(j-1)]) + 1

    }
    
    k <- sum(ns[1:j])

    mn <- apply(x[i:k,], 2, mean)

    s <- do.call(cov_fun, list(x[i:k, ]))

    rbind(mn, s)
  }

  u <- do.call(cbind, result)

  mn.cr <- combine_cov_estimates(u, ns)

  result <- mn.cr[2:NROW(mn.cr), ]

  rownames(result) <- colnames(x)

  colnames(result) <- colnames(x)

  result
}