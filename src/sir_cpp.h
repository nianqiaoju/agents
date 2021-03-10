#ifndef _INCL_SIR_
#define _INCL_SIR_
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;

// this is a header file for all functions releated to the sir model

NumericVector sir_bif_cpp(const IntegerVector & y,
                          const IntegerMatrix & nexts,
                          const IntegerMatrix & nexti,
                          const NumericMatrix & logfbar,
                          const double & rho,
                          const int & N);


List create_fbar_matrix(const double lambda, 
                        const double gamma,
                        const int N,
                        const double c);

void update_fbar_matrix(IntegerMatrix & nexts,
                        IntegerMatrix & nexti, 
                        NumericMatrix & logfbar, 
                        const double & lambda,
                        const double & gamma, 
                        const int & N,
                        const double & c);



int sir_get_state_index(const int s, const int i, const int N){
  int index = (N - s) * (N- s + 1) / 2 + N - s - i;
  return(index);
}

double sir_logfbar(const int snow, const int inow, const int snext, const int inext,
               const double lambda, const double gamma, const int N){
  return(R::dbinom(snow + inow - inext - snext, inow, gamma, true) + R::dbinom(snow - snext, snow, lambda * inow / N, true));
};

#endif