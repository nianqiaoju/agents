#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sis_get_alpha_full_cpp(const LogicalVector &xx,
                                const NumericVector &lambda, 
                                const NumericVector &gamma){
  NumericVector alpha (xx.length());
  double num_infections = sum(xx);
  // Rcout << "there are " << num_infections << "infections\n";
  double prop_infections = num_infections / xx.length();
  for(int n = 0; n < xx.length(); n++){
    alpha[n] = (xx[n]) ? (1 - gamma[n]) : (prop_infections * lambda[n]);
  }
  return alpha;
}
