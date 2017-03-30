#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector distributed_var_base(NumericVector v1, NumericVector v2, int n_A, int n_B) {  
  double delta = v2[0] - v1[0];
  double m_a = v1[1] * (n_A - 1);
  double m_b = v2[1] * (n_B - 1);
  double M2 = m_a + m_b + delta * delta * n_A * n_B / (n_A + n_B);
  NumericVector res;
  res.push_back((v1[0]*n_A + v2[0]*n_B)/(n_A + n_B));
  res.push_back(M2/(n_A + n_B - 1));
  return res;
}

// [[Rcpp::export]]
NumericVector custom_mean_var(NumericVector x) { 
  NumericVector y;

  int n = 0;
  double mn = 0;
  double M2 = 0;
  
  for(int i = 0; i < x.size() ; i++) {
    double xi = x[i];
    n++;
    double delta = xi - mn;
    mn += delta/n;
    double delta2 = xi - mn;
    M2 += delta * delta2;
  }
  
  if (n < 2) {
    y.push_back(x[0]);
    y.push_back(0);
  } else {
    y.push_back(mn);
    y.push_back(M2/(n-1));
  }
  
  return y;
}

// [[Rcpp::export]]
NumericVector combine_stats(NumericVector mns, NumericVector vars, NumericVector ns) {
  if (mns.size() == 2) {
    NumericVector res;
    NumericVector v1;
    NumericVector v2;
    v1.push_back(mns[0]);
    v1.push_back(vars[0]);
    v2.push_back(mns[1]);
    v2.push_back(vars[1]);
    
    return(distributed_var_base(v1, v2, ns[0], ns[1]));
  }
  
  NumericVector mns_left;
  for (int i = 0; i < mns.size()/2 ; i++) {
    mns_left.push_back(mns[i]);
  }
  NumericVector mns_right;
  for (int i = mns.size()/2; i < mns.size() ; i++) {
    mns_right.push_back(mns[i]);
  }

  NumericVector vars_left;
  for (int i = 0; i < vars.size()/2 ; i++) {
    vars_left.push_back(vars[i]);
  }
  NumericVector vars_right;
  for (int i = vars.size()/2; i < vars.size() ; i++) {
    vars_right.push_back(vars[i]);
  }

  NumericVector ns_left;
  for (int i = 0; i < ns.size()/2 ; i++) {
    ns_left.push_back(ns[i]);
  }
  NumericVector ns_right;
  for (int i = ns.size()/2; i < ns.size() ; i++) {
    ns_right.push_back(ns[i]);
  }  
  
  NumericVector v1 = combine_stats(mns_left, vars_left, ns_left);
  NumericVector v2 = combine_stats(mns_right, vars_right, ns_right);
  
  int ns_sum_left = 0;
  for (int i = 0; i < ns_left.size() ; i++) {
    ns_sum_left += ns_left[i];
  }
  int ns_sum_right = 0;
  for (int i = 0; i < ns_right.size() ; i++) {
    ns_sum_right += ns_right[i];
  }
  
  NumericVector res = distributed_var_base(v1, v2, ns_sum_left, ns_sum_right);
  return res;
}