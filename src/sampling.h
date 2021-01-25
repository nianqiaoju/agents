#ifndef _INCL_SAMPLING_
#define _INCL_SAMPLING_
#include <RcppEigen.h>
#include <string>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;

// this is a header file for all functions releated sampling or resampling
NumericVector get_runifs(const int n);

LogicalVector idchecking_cpp(const int sum_x,
                   const NumericVector &alpha,
                   const NumericVector &random_number);

// sample from multinomial distribution given log weights and a uniform draw
int multinomial_cpp(const NumericVector &logweights, double uniform);

IntegerVector multinomial_resampling_cpp(const NumericVector & weights, int ndraws, const NumericVector & rand);

#endif