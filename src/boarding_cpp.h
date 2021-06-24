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

int boarding_lowdim2index(const int N, const int scnt, const int icnt);

void boarding_sample_x_given_si_sparse(IntegerMatrix & xx,
                                    const NumericMatrix & alphas2i,
                                    const NumericMatrix & alphai2i, 
                                    const NumericVector & lambda,
                                    const NumericVector & gamma, 
                                    const IntegerVector & snext,
                                    const IntegerVector & inext,
                                    const int & N,
                                    const int & P);

#endif