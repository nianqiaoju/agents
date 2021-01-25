#include <Rcpp.h>
#include "logdensities.h"
#include "sis_cpp.h"
using namespace Rcpp;
using namespace std;

// this function performs a complete single-site scan for all t and all n.

// [[Rcpp::export]]


void sis_xx_gibbs_singlesite_full_cpp(LogicalMatrix xx, 
                                      const IntegerVector &y,
                                      const NumericVector &alpha0,
                                      const NumericVector &lambda, 
                                      const NumericVector &gamma, 
                                      const double rho){
  GetRNGstate();
  
  int population_size = xx.ncol();
  LogicalVector xxt(population_size);
  NumericVector alphaprev2t(population_size);
  NumericVector alphat2next(population_size);
  LogicalVector xxnext(population_size);
  double lw0, lw1, p1, maxlw;
  
  int time_step = 0;
  // Rcout << "t = " << time_step << "\n";
  int infection_t = 0;
  alphaprev2t = alpha0;
  xxt = xx(_ ,time_step);
  xxnext = xx(_, time_step + 1);
  for(int n = 0; n < population_size; n++){
    xxt[n] = false;
    infection_t = sum(xxt);
    alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
    lw0 = log(1- alphaprev2t[n]) + R::dbinom(y[time_step], infection_t, rho, true) + logdbern_sum_cpp(alphat2next,xxnext);
    xxt[n] = true;
    alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
    lw1 = log(alphaprev2t[n]) + R::dbinom(y[time_step], infection_t + 1, rho, true) +  logdbern_sum_cpp(alphat2next, xxnext);
    // normalize probabilities
    maxlw = max(lw0,lw1);
    p1 = exp(lw1- maxlw) / (exp(lw1 - maxlw) + exp(lw0 - maxlw));
    // sample xt[n]
    if(unif_rand() < p1){
      xxt[n] = true;
    }else{
      xxt[n] = false;
    }
  }
  xx(_, time_step) = xxt;
  for(time_step = 1; (time_step + 1) < y.length(); time_step++){
    xxt = xx(_, time_step);
    alphaprev2t = sis_get_alpha_full_cpp(xx(_, time_step - 1), lambda, gamma);
    xxnext = xx(_, time_step + 1);
    // Rcout << "BEFORE: t = " << time_step << "y = " <<  y[time_step] << "infection count is " <<  sum(xxt) << "\n";
    for(int n = 0; n < population_size; n++){
      xxt[n] = false;
      infection_t = sum(xxt);
      alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
      lw0 = log(1 - alphaprev2t[n]) + R::dbinom(y[time_step], infection_t, rho, true) + logdbern_sum_cpp(alphat2next, xxnext);
      xxt[n] = true;
      alphat2next = sis_get_alpha_full_cpp(xxt,lambda,gamma);
      lw1 = log(alphaprev2t[n]) + R::dbinom(y[time_step], infection_t + 1, rho, true) + logdbern_sum_cpp(alphat2next, xxnext);
      // normalize probabilities
      maxlw = max(lw0,lw1);
      p1 = exp(lw1- maxlw) / (exp(lw1 - maxlw) + exp(lw0 - maxlw));
      // sample xt[n]
      if(unif_rand() < p1){
        xxt[n] = true;
      }else{
        xxt[n] = false;
      }
    }
    // Rcout << "AFTER: t = " << time_step << "y = " <<  y[time_step] << "infection count is " <<  sum(xxt) << "\n";
    xx(_, time_step) = xxt;
  }
  // final step
  // Rcout << "t = " << time_step << "\n";
  xxt = xx(_, time_step);
  alphaprev2t = sis_get_alpha_full_cpp(xx(_, time_step - 1), lambda, gamma);
  for(int n = 0; n < population_size; n++){
    xxt[n] = false;
    infection_t = sum(xxt);
    lw0 = log(1 - alphaprev2t[n]) + R::dbinom(y[time_step], infection_t, rho, true);
    xxt[n] = true;
    lw1 = log(alphaprev2t[n]) + R::dbinom(y[time_step], infection_t + 1, rho, true);
    // normalize probabilities
    maxlw = max(lw0,lw1);
    p1 = exp(lw1- maxlw) / (exp(lw1 - maxlw) + exp(lw0 - maxlw));
    // sample xt[n]
    if(unif_rand() < p1){
      xxt[n] = true;
    }else{
      xxt[n] = false;
    }
  }
  xx(_, time_step) = xxt;
  PutRNGstate();
}
