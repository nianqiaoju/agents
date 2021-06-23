#ifndef _INCL_BOARDING_
#define _INCL_BOARDING_
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

/*
 * Functions in this header file are written for the boarding school example with N = 120..
 */
NumericVector smallpox_bif_create_cpp(const IntegerVector & y,
                                      const IntegerMatrix & nexts,
                                      const IntegerMatrix & nexti,
                                      const NumericMatrix & logfbar,
                                      const double & rho,
                                      const int & N);

void smallpox_bif_update_cpp(NumericVector & logpolicy, 
                             IntegerMatrix & nexts,
                             IntegerMatrix & nexti,
                             NumericMatrix & logfbar,
                             const double & lambda,
                             const double & gamma, 
                             const double & rho,
                             const IntegerVector & y,
                             const int & N,
                             const double & c);


#endif