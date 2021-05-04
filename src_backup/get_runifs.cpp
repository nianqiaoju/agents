#include <Rcpp.h>
#include "sampling.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_runifs(const int n){
  GetRNGstate();
  NumericVector us (n);
  for(int i = 0; i < n; i++){
    us[i] = unif_rand();
  }
  PutRNGstate();
  return us;
}