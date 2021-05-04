#ifndef _INCL_SIS_
#define _INCL_SIS_
#include <RcppEigen.h>
#include "logsumexp.h"
#include "logdensities.h"
#include "sampling.h"
#include <string>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;

// this is a header file for all functions releated to the sis model

NumericMatrix sis_forward_algorithm_cpp(NumericVector logf0, 
                                        NumericMatrix logdtransition,  
                                        NumericVector y, 
                                        NumericVector all_sum_x, 
                                        double rho);

// the following functions computes alpha given parameters and current state, assuming a fully connected network

// compute alpha(xx) given agent states and parameters
NumericVector sis_get_alpha_full_cpp(const LogicalVector &xx,
                                     const NumericVector &lambda, 
                                     const NumericVector &gamma);

// gibbs samplers
void sis_xx_gibbs_swap_full_cpp(LogicalMatrix xx, 
                                const NumericVector &alpha0,
                                const NumericVector &lambda, 
                                const NumericVector &gamma);

// this function modifies the LogicalMatrix xx directly,hence returns void.
void sis_xx_gibbs_singlesite_full_cpp(LogicalMatrix xx, 
                                      const IntegerVector &y,
                                      const NumericVector alpha0,
                                      const NumericVector lambda, 
                                      const NumericVector gamma, 
                                      const double rho);

// blocked gibbs updates
void sis_xx_gibbs_blocked_full_cpp(LogicalMatrix  xx,
                                   const IntegerVector &y,
                                   const NumericVector &alpha0,
                                   const NumericVector &lambda,
                                   const NumericVector &gamma,
                                   const double rho,
                                   const int block_size,
                                   const LogicalMatrix &state_space);

// // particle filters
// returns estimate of log marginal likelihood
// implementations are specific to fully-connected networks

// double sis_bpf_full_cpp(const IntegerVector &y,
//   const NumericVector &alpha0,
//   const NumericVector &lambda,
//   const NumericVector &gamma,
//   const double rho,
//   const int num_particles,
//   const double threshold,
//   const int population_size
//   );

double sis_apf_exact_full_cpp(const IntegerVector &y,
  const NumericVector &alpha0,
  const NumericVector &lambda,
  const NumericVector &gamma,
  const double rho,
  const int num_particles,
  const double threshold,
  const int population_size
  );

// NumericMatrix sis_bif_cpp_tp(const IntegerVector &y,
//   const double lambda_bar,
//   const double gamma_bar);
// 
// double sis_csmc_full_cpp(const IntegerVector &y,
//   const NumericVector &alpha0,
//   const NumericVector &lambda,
//   const NumericVector &gamma,
//   const double rho,
//   const int num_particles,
//   const double threshold,
//   const int population_size,
//   const NumericMatrix policy
//   );

#endif