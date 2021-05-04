#include <Rcpp.h>
#include "sampling.h"
using namespace Rcpp;

// given weights (normalized or unnormalied)
// sample without replacement ndraws
// given random numbers
// output is zero based

// [[Rcpp::export]]
IntegerVector multinomial_resampling_cpp(const NumericVector & weights, int ndraws, const NumericVector & rand){
  IntegerVector ancestors(ndraws);
  NumericVector cumsumw = cumsum(weights);
  int i;
  for (int n = 0; n < ndraws; n++){
    i = 0;
    while (cumsumw(i) < rand(n)){
      i = i + 1;
    }
    ancestors(n) = i;
  }
  return ancestors;
}