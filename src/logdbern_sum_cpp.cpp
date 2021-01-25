#include <RcppEigen.h>
#include <Rcpp.h>
#include "logdensities.h"
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;


// assume that x is a binary vector such that each entry is independently Bernoulli(alpha[n])
// compute the log joint probability mass function of an observation x

// [[Rcpp::export]]
double logdbern_sum_cpp (const NumericVector & alpha, const  LogicalVector &x){
  // if alpha is zero vector and x is zero vector, then probability is one
  double lprod = 0;
  int N = x.size();
  for (int n = 0; n < N; n++){
    lprod += ((x[n] == 1)? log(alpha[n]) : log(1 - alpha[n]));
  }
  return(lprod);
}

