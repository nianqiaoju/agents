#include <Rcpp.h>
#include "smallpox.h"
#include "sir_cpp.h"
#include "logdensities.h"
#include "sampling.h"
using namespace Rcpp;
using namespace std;

/*
 * create the logpolicy matrix for the smallpox SIR model
 * where observations are Binomial(Rt, rho)
 * this function should run after nexts, nexti and logfbar are created in sir_creat_f_matrix function
 */

// [[Rcpp::export]]
NumericVector smallpox_bif_create_cpp(const IntegerVector & y, 
                                      const IntegerMatrix & nexts,
                                      const IntegerMatrix & nexti,
                                      const NumericMatrix & logfbar,
                                      const double & rho, 
                                      const int & N) {
  int days = y.size();
  // Rcout << "observation length = " <<  days <<"\n";
  // initialize a vector to -Inf
  NumericVector logpolicy((N+1) * (N+1)  * days);
  std::fill(logpolicy.begin(), logpolicy.end(), R_NegInf);
  // deal with the terminal time first
  int t = days - 1;
  // Rprintf("at day t = %i\n", t);
  for(int inow = N - y[t]; inow >= 0; inow--){
    for(int snow = 0; snow + inow <= N - y[t]; snow++){
      // Rprintf("s = %i, i = %i\n", snow, inow);
      logpolicy[t * (N+1) * (N+1) + (N+1) * inow + snow] = R::dbinom(y[t], N - inow - snow, rho, true);
    }
  }
  NumericVector logg(N + 1);
  double curr_sum, add_item, curr_max;
  int cntstatenext;
  int istatenow;
  for(int t = days - 2; t >= 0; t--){
    // Rprintf("at day t = %i\n", t);
    // update observation densities at time t
    // yt ~ Binom(rt, rho)
    std::fill(logg.begin(), logg.end(), R_NegInf);
    for(int r = y[t]; r <= N; r++){
      logg[r] = R::dbinom(y[t], r, rho, true);
    }
    
    // compute the conditional expecation for each value of st and it
    for(int snow = N - y[t]; snow >= 0; snow--){
      istatenow = sir_get_state_index(snow, N - snow - y[t], N); 
      for(int inow = N - snow - y[t]; inow >= 0; inow--){
        curr_sum = R_NegInf; // initialize the conditonal expectation to 0
        // can also set curr_sum to a very small value to satisfy the sufficient support condition
        cntstatenext = 0;
        while(cntstatenext < logfbar.ncol() && nexti(istatenow, cntstatenext) >= 0){
          add_item = logfbar(istatenow, cntstatenext) + logpolicy[(t + 1) * (N+1) * (N+1) + (N+1) * nexti(istatenow, cntstatenext) + nexts(istatenow, cntstatenext)];
          if(!Rcpp::traits::is_infinite<REALSXP>(add_item)){
            curr_max = max(curr_sum, add_item);
            curr_sum = log(exp(curr_sum - curr_max) + exp(add_item - curr_max)) + curr_max;
          }
          cntstatenext++;
        }
        // psi_t(st,it) is observation density * conditional expectation
        logpolicy[t * (N+1) * (N+1) + (N+1) * (inow) + snow] = curr_sum + logg[N - inow - snow];
        istatenow++;
      }
    }
  }
  // convert vector to  3d-array
  IntegerVector dim = {N+1, N+1, days};
  logpolicy.attr("dim") = dim;
  return logpolicy;
}


/* For the smallpox dataset, 
 * update the backward information filter given new values of lambda, gamma and rho.
 * this function should be run after smallpox_bif_create_cpp and create_fbar_matrix, 
 * and it modifies logpolicy and logfbar matrix directly.
 */

// [[Rcpp::export]]
void smallpox_bif_update_cpp(NumericVector & logpolicy, 
                             IntegerMatrix & nexts,
                             IntegerMatrix & nexti, 
                             NumericMatrix & logfbar, 
                             const double & lambda,
                             const double & gamma, 
                             const double & rho,
                             const IntegerVector & y, 
                             const int & N,
                             const double & c){
  update_fbar_matrix(nexts, nexti, logfbar, lambda, gamma, N, c);
  std::fill(logpolicy.begin(), logpolicy.end(), R_NegInf);
  // start at the terminal time first
  int days = y.size();
  int t = days - 1;
  for(int inow = N - y[t]; inow >= 0; inow--){
    for(int snow = 0; snow + inow <= N - y[t]; snow++){
      // Rprintf("s = %i, i = %i\n", snow, inow);
      logpolicy[t * (N+1) * (N+1) + (N+1) * inow + snow] = R::dbinom(y[t], N - inow - snow, rho, true);
    }
  }
  
  NumericVector logg(N + 1);
  double curr_sum, add_item, curr_max;
  int cntstatenext;
  int istatenow;
  for(int t = days - 2; t >= 0; t--){
    // update observation densities at time t
    // yt ~ Binom(rt, rho)
    std::fill(logg.begin(), logg.end(), R_NegInf);
    for(int r = y[t]; r <= N; r++){
      logg[r] = R::dbinom(y[t], r, rho, true);
    }
    // compute the conditional expecation for each value of st and it
    for(int snow = N - y[t]; snow >= 0; snow--){
      istatenow = sir_get_state_index(snow, N - snow - y[t], N); 
      for(int inow = N - snow - y[t]; inow >= 0; inow--){
        curr_sum = R_NegInf; // initialize the conditonal expectation to 0
        // can also set curr_sum to a very small value to satisfy the sufficient support condition
        cntstatenext = 0;
        while(cntstatenext < logfbar.ncol() && nexti(istatenow, cntstatenext) >= 0){
          add_item = logfbar(istatenow, cntstatenext) + logpolicy[(t + 1) * (N+1) * (N+1) + (N+1) * nexti(istatenow, cntstatenext) + nexts(istatenow, cntstatenext)];
          if(!Rcpp::traits::is_infinite<REALSXP>(add_item)){
            curr_max = max(curr_sum, add_item);
            curr_sum = log(exp(curr_sum - curr_max) + exp(add_item - curr_max)) + curr_max;
          }
          cntstatenext++;
        }
        logpolicy[t * (N+1) * (N+1) + (N+1) * (inow) + snow] = curr_sum + logg[N - inow - snow];
        istatenow++;
      }
    }
  }
}
