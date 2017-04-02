#include <Rcpp.h>

using namespace Rcpp;

double custom_covar(double* x1, double* x2, int size) { 
  int n = 0;
  double mn1 = 0;
  double mn2 = 0;
  double M12 = 0;
  
  for(int i = 0; i < size ; i++) {
    double xi1 = x1[i];
    double xi2 = x2[i];
    n++;
    double delta1 = (xi1 - mn1)/n;
    mn1 += delta1;
    double delta2 = (xi2 - mn2)/n;
    mn2 += delta2;
    M12 += (n - 1) * delta1 * delta2 - M12/n;
  }
  
  if (n < 2) {
    return 0;
  } else {
    return (M12 *n/(n-1));
  }
}

double custom_covar(double* x1, double* x2, double mn1, double mn2, int size) { 
  int n = size;
  double M12 = 0;
  
  for(int i = 0; i < size ; i++) {
    double delta1 = (x1[i] - mn1);
    double delta2 = (x2[i] - mn2);
    M12 += (delta1 * delta2)/(n-1);
  }
  
  if (n < 2) {
    return 0;
  } else {
    return (M12);
  }
}

//' custom_covar
//' 
//' @param x1 ...
//' @param x2 ...
//' 
//' @export
//' 
// [[Rcpp::export]]
double custom_covar(NumericVector x1, NumericVector x2) { 
  return(custom_covar(&x1[0], &x2[0], x1.size()));
  
}

//' combine_covs
//' 
//' @param mn_covs ...
//' @param ns ...
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericMatrix custom_multi_covar(NumericMatrix s) { 
  NumericMatrix res(s.ncol(), s.ncol());
  double* Means = new double[s.ncol()];
  double* sx = &s[0];
  for (int i = 0 ; i < s.ncol() ; i++) {
    double sm = 0;
    for (int j = 0 ; j < s.nrow() ; j++) {
      sm += s(j, i);
    }
    sm /= s.nrow();
    Means[i] = sm;
  }
  
  int n = s.nrow();
  for (int i = 0 ; i < s.ncol() ; i++) {
    for (int j = 0 ; j <= i ; j++) {
      double M12 = 0;
      double mn1 = Means[i];
      double mn2 = Means[j];
      for(int k = 0; k < n ; k++) {
        M12 += (sx[i*n + k] - mn1) * (sx[j*n + k] - mn2)/(n-1);   // division inside
                                                                  // the loop makes it slower
                                                                  // but numerically stabler
      }
      res(i, j) = M12;
      res(j, i) = res(i, j);
    }
  }
  delete[] Means;
  return res;
}

//' combine_covs_base
//' 
//' @param mn_covs1 ...
//' @param mn_covs2 ...
//' @param ns1 ...
//' @param ns2 ...
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericMatrix combine_covs_base(NumericMatrix mn_covs1, NumericMatrix mn_covs2, int ns1, int ns2) {
  NumericMatrix res(mn_covs1.nrow(), mn_covs1.nrow()-1);
  for (int i = 1 ; i < mn_covs1.nrow() ; i++) {
    res(0, i - 1) = (mn_covs1(0, i - 1) * ns1 + mn_covs2(0, i - 1) * ns2)/(ns1 + ns2);
    for (int j = 0 ; j < mn_covs1.nrow() - 1 ; j++) {
      double cx = ((ns1-1) * mn_covs1(i, j) + (ns2-1) * mn_covs2(i, j) + (mn_covs1(0, i-1) - mn_covs2(0, i - 1)) * (mn_covs1(0, j) - mn_covs2(0, j)) * (ns1) * (ns2)/(ns1 + ns2))/(ns1 + ns2- 1);
      res(i, j) = cx;
    }
  }
  
  return res;
}

//' combine_covs
//' 
//' @param mn_covs ...
//' @param ns ...
//' 
//' @export
//' 
// [[Rcpp::export]]
NumericMatrix combine_covs(NumericMatrix mn_covs, NumericVector ns) {
  if (ns.size() == 2) {
    return combine_covs_base(mn_covs(_,Range(0, mn_covs.ncol()/2 - 1)), mn_covs(_,Range(mn_covs.ncol()/2, mn_covs.ncol()-1)), ns[0], ns[1]);
  } else {
    NumericMatrix c_left = mn_covs(_,Range(0, mn_covs.ncol()/2-1));
    NumericMatrix c_right = mn_covs(_,Range(mn_covs.ncol()/2, mn_covs.ncol()-1));
    NumericVector ns_left(ns.size()/2);
    for (int i = 0 ; i < ns.size()/2 ; i++) {
      ns_left[i] = ns[i];
    }
    
    NumericVector ns_right(ns.size()/2);
    for (int i = ns.size()/2 ; i < ns.size() ; i++) {
      ns_right[i - ns.size()/2] = ns[i];
    }
    
    NumericMatrix res_left = combine_covs(c_left, ns_left);
    NumericMatrix res_right = combine_covs(c_right, ns_right);
    int ns_left_sum = 0;
    for (int i = 0 ; i < ns_left.size() ; i++) {
      ns_left_sum += ns_left[i];
    }
    int ns_right_sum = 0;
    for (int i = 0 ; i < ns_right.size() ; i++) {
      ns_right_sum += ns_right[i];
    }
    
    return(combine_covs_base(res_left, res_right, ns_left_sum, ns_right_sum));
  }
}