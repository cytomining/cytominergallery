// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// custom_covar
double custom_covar(NumericVector x1, NumericVector x2);
RcppExport SEXP cytominergallery_custom_covar(SEXP x1SEXP, SEXP x2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x2(x2SEXP);
    rcpp_result_gen = Rcpp::wrap(custom_covar(x1, x2));
    return rcpp_result_gen;
END_RCPP
}
// custom_multi_covar
NumericMatrix custom_multi_covar(NumericMatrix s);
RcppExport SEXP cytominergallery_custom_multi_covar(SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(custom_multi_covar(s));
    return rcpp_result_gen;
END_RCPP
}
// combine_covs_base
NumericMatrix combine_covs_base(NumericMatrix mn_covs1, NumericMatrix mn_covs2, int ns1, int ns2);
RcppExport SEXP cytominergallery_combine_covs_base(SEXP mn_covs1SEXP, SEXP mn_covs2SEXP, SEXP ns1SEXP, SEXP ns2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mn_covs1(mn_covs1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mn_covs2(mn_covs2SEXP);
    Rcpp::traits::input_parameter< int >::type ns1(ns1SEXP);
    Rcpp::traits::input_parameter< int >::type ns2(ns2SEXP);
    rcpp_result_gen = Rcpp::wrap(combine_covs_base(mn_covs1, mn_covs2, ns1, ns2));
    return rcpp_result_gen;
END_RCPP
}
// combine_covs
NumericMatrix combine_covs(NumericMatrix mn_covs, NumericVector ns);
RcppExport SEXP cytominergallery_combine_covs(SEXP mn_covsSEXP, SEXP nsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mn_covs(mn_covsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ns(nsSEXP);
    rcpp_result_gen = Rcpp::wrap(combine_covs(mn_covs, ns));
    return rcpp_result_gen;
END_RCPP
}