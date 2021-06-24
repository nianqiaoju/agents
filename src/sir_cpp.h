#ifndef _INCL_SIR_
#define _INCL_SIR_
#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;

// this is a header file for all functions releated to the sir model

NumericVector sir_bif_create(const IntegerVector & y,
                          const IntegerMatrix & nexts,
                          const IntegerMatrix & nexti,
                          const NumericMatrix & logfbar,
                          const double & rho,
                          const int & N);

void sir_bif_update(NumericVector & logpolicy, 
                        IntegerMatrix & nexts,
                        IntegerMatrix & nexti,
                        NumericMatrix & logfbar,
                        const double & lambda,
                        const double & gamma, 
                        const double & rho,
                        const IntegerVector & y,
                        const int & N,
                        const double & c);


List sir_fbar_create(const double lambda, 
                        const double gamma,
                        const int N,
                        const double c);

void sir_fbar_update(IntegerMatrix & nexts,
                        IntegerMatrix & nexti, 
                        NumericMatrix & logfbar, 
                        const double & lambda,
                        const double & gamma, 
                        const int & N,
                        const double & c);

int sir_get_state_index(const int s, const int i, const int N);


double sir_logfbar(const int & snow,
                   const int & inow, 
                   const int & snext, 
                   const int & inext,
                   const double & lambda, 
                   const double & gamma, 
                   const int & N);

void sir_logf_update(NumericMatrix & logf,
                         const IntegerVector & xxprev,
                         const NumericVector & lambda_v,
                         const NumericVector & gamma_v,
                         const int & N);


void sir_alpha_update_sparse(NumericMatrix alphas2i, 
                             NumericMatrix alphai2i, 
                             const IntegerMatrix & xts,
                             const NumericVector & lambda,
                             const NumericVector  & gamma, 
                             const IntegerMatrix  & neighbors,
                             const int & N);
  


IntegerMatrix sir_sample_x_given_si(IntegerMatrix & xts,
                           const NumericVector & lambda,
                           const NumericVector & gamma, 
                           const IntegerVector & scount,
                           const IntegerVector & icount,
                           const int & N,
                           const int & P);
bool xory(bool x, bool y);


#endif
