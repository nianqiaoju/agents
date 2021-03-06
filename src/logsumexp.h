#ifndef _INCL_LOGSUMEXP_
#define _INCL_LOGSUMEXP_
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;

// takes a vector of log weights, and compute the log of the sum of the weights
double lw_logsum_cpp (NumericVector lw);

// takes a vectors of log weights and normalize the weights
NumericVector lw_normalize_cpp(const NumericVector & lw);

#endif