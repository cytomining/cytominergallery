// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// online_covar
double online_covar(NumericVector x1, NumericVector x2);
RcppExport SEXP cytominergallery_online_covar(SEXP x1SEXP, SEXP x2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x2(x2SEXP);
    rcpp_result_gen = Rcpp::wrap(online_covar(x1, x2));
    return rcpp_result_gen;
END_RCPP
}
// two_pass_multi_covar
NumericMatrix two_pass_multi_covar(NumericMatrix s);
RcppExport SEXP cytominergallery_two_pass_multi_covar(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(two_pass_multi_covar(s));
    return rcpp_result_gen;
END_RCPP
}
// combine_cov_estimates
NumericMatrix combine_cov_estimates(NumericMatrix batch_mean_cov, NumericVector b);
RcppExport SEXP cytominergallery_combine_cov_estimates(SEXP batch_mean_covSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type batch_mean_cov(batch_mean_covSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(combine_cov_estimates(batch_mean_cov, b));
    return rcpp_result_gen;
END_RCPP
}
