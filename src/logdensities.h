#ifndef _INCL_LOG_DENSITIES_
#define _INCL_LOG_DENSITIES_
#include <RcppEigen.h>
#include "logsumexp.h"
#include <string>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;

// this is a header file for all functions releated to log densities

double logdbern_sum_cpp (const NumericVector &alpha, const LogicalVector &x);

NumericVector logdpoisbinom_cpp(const NumericVector &alpha);

#endif