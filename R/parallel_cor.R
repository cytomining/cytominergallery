#' Parallel correlation function
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
parallel_cor <- function(x, splits = 2, cores = 2, cov_fun = "two_pass_multi_covar") {
  
  cov_mat <- parallel_cov(x, splits = splits, cores = cores, cov_fun = cov_fun)
  
  nrm <- diag(diag(cov_mat)^-0.5)
  
  cov_mat <- cov_mat %>% as.matrix()
  
  nrm <- as.matrix(nrm)
  
  result <- nrm %*% cov_mat %*% nrm
  
  rownames(result) <- colnames(x)
  
  colnames(result) <- colnames(x)
  
  result
}