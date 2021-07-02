#ifndef _INCL_BOARDING_
#define _INCL_BOARDING_
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

/*
 * Functions in this header file are written for the boarding school example with N = 763.
 */

IntegerMatrix boarding_all_lowdim_states(const int N);

int boarding_lowdim2index(const int N, const int scnt, const int icnt);

NumericMatrix boarding_bif_create(const IntegerVector & y,
                                  const IntegerMatrix & all_lowdim_states,
                                  const double & lambda, 
                                  const double & gamma, 
                                  const double & rho,
                                  const int & N,
                                  const double & c);

void boarding_bif_update(NumericMatrix logpolicy, 
                         const IntegerVector & y, 
                         const IntegerMatrix & all_lowdim_states,
                         const double & lambda, 
                         const double & gamma, 
                         const double & rho,
                         const int & N,
                         const double &c);


List boarding_fbar_create(const double lambda, 
                          const double gamma,
                          const int N,
                          const double c);

void boarding_fbar_update(NumericMatrix logfbar,
                          IntegerMatrix nextsi,
                          const double lambda,
                          const double gamma,
                          const double N,
                          const double c);

NumericMatrix boarding_bif_create_fast(const IntegerVector & y,
                                       const IntegerMatrix & nextsi,
                                       const NumericMatrix & logfbar,
                                       const double & rho,
                                       const int & N,
                                       const double & c);

void boarding_bif_update_fast(NumericMatrix logpolicy,
                              const IntegerVector & y,
                              const NumericMatrix & logfbar,
                              const IntegerMatrix & nextsi,
                              const double & rho,
                              const int & N,
                              const double & c);


void boarding_logf_update_sparse(NumericMatrix logf,
                                 NumericMatrix alphas2i,
                                 NumericMatrix alphai2i,
                                 const IntegerMatrix  & xts,
                                 const double & lambda,
                                 const double & gamma,
                                 const IntegerMatrix & neighbors,
                                 const int & N);


void boarding_logf_update_full(NumericMatrix logf,
                                 NumericMatrix alphas2i,
                                 NumericMatrix alphai2i,
                                 const IntegerMatrix  & xts,
                                 const double & lambda,
                                 const double & gamma,
                                 const int & N);



void boarding_sample_x_given_si(IntegerMatrix & xts,
                                const NumericMatrix & alphas2i,
                                const NumericMatrix & alphai2i, 
                                const IntegerVector & snext,
                                const IntegerVector & inext,
                                const int & N,
                                const int & P);

#endif