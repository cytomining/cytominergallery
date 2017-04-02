#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector custom_mean_var(NumericVector x) { 
  NumericVector y;
  double* x1 = &x[0];
  int n = 0;
  double mn = 0;
  double M2 = 0;
  
  for(int i = 0; i < x.size() ; i++) {
    double xi = x1[i];
    n++;
    double delta = xi - mn;
    mn += delta/n;
    double delta2 = xi - mn;
    M2 += delta * delta2;
  }
  
  if (n < 2) {
    y.push_back(x1[0]);
    y.push_back(0);
  } else {
    y.push_back(mn);
    y.push_back(M2/(n-1));
  }
  
  return y;
}

double* distributed_var_base(double mn1, double mn2, double var1, double var2, int n_A, int n_B) {  
  double delta = mn1 - mn2;
  double m_a = var1 * (n_A - 1);
  double m_b = var2 * (n_B - 1);
  double M2 = m_a + m_b + delta * delta * n_A * n_B / (n_A + n_B);
  double* res = new double[2];
  res[0] = ((mn1*n_A + mn2*n_B)/(n_A + n_B));
  res[1] = (M2/(n_A + n_B - 1));
  return res;
}

double* combine_stats_var_internal(double* mns, double* vars, double* ns, int mns_len, int vars_len, int ns_len) {
  if (mns_len == 2) {
    return(distributed_var_base(mns[0], mns[1], vars[0], vars[1], ns[0], ns[1]));
  }
  
  double* v1 = combine_stats_var_internal(mns, vars, ns, mns_len/2, vars_len/2, ns_len/2);
  double* v2 = combine_stats_var_internal(mns + mns_len/2, vars + vars_len/2, ns + ns_len/2, mns_len/2, vars_len/2, ns_len/2);
  
  int ns_sum_left = 0;
  for (int i = 0; i < ns_len/2 ; i++) {
    ns_sum_left += ns[i];
  }
  int ns_sum_right = 0;
  for (int i = ns_len/2; i < ns_len ; i++) {
    ns_sum_right += ns[i];
  }
  
  double* res = distributed_var_base(v1[0], v2[0], v1[1], v2[1], ns_sum_left, ns_sum_right);
  return res;
}

// [[Rcpp::export]]
NumericVector combine_stats_var(NumericVector mns, NumericVector vars, NumericVector ns) {
  double* d = (combine_stats_var_internal(&mns[0], &vars[0], &ns[0], mns.size(), vars.size(), ns.size()));
  NumericVector res;
  res.push_back(d[0]);
  res.push_back(d[1]);
  return res;
}