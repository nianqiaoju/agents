#ifndef _INCL_BOARDING_
#define _INCL_BOARDING_
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

/*
 * Functions in this header file are written for the boarding school example with N = 763.
 */

NumericMatrix boarding_bif_create_cpp(const IntegerVector & y,
                                      const double & lambda, 
                                      const double & gamma, 
                                      const double & rho,
                                      const int & N,
                                      const double & c); // current implementation assumes network_type == "full"

#endif