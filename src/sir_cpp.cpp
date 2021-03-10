#include <Rcpp.h>
#include "sir_cpp.h"
using namespace Rcpp;
using namespace std;

// perform the BIF approximation for the SIR model
// after the approximate trasition kernel fbar has been created by create_fbar_matrix 
// or modified in update_fbar_matrix

// [[Rcpp::export]]
NumericVector sir_bif_cpp(const IntegerVector & y, 
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
  for(int inow = N; inow >= y[t]; inow--){
    for(int snow = 0; snow + inow <= N; snow++){
      // Rprintf("s = %i, i = %i\n", snow, inow);
      logpolicy[t * (N+1) * (N+1) + (N+1) * inow + snow] = R::dbinom(y[t], inow, rho, true);
    }
  }
  NumericVector logg(N + 1);
  double curr_sum, add_item, curr_max;
  int cntstatenext;
  int istatenow;
  for(int t = days - 2; t >= 0; t--){
    // Rprintf("at day t = %i\n", t);
    // update observation densities at time t
    std::fill(logg.begin(), logg.end(), R_NegInf);
    for(int i = y[t]; i <= N; i++){
      logg[i] = R::dbinom(y[t], i, rho, true);
    }
    
    // compute the conditional expecation for each value of st and it
    for(int snow = N - y[t]; snow >= 0; snow--){
      istatenow = sir_get_state_index(snow, N - snow, N); 
      for(int inow = N - snow; inow >= y[t]; inow--){
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
        logpolicy[t * (N+1) * (N+1) + (N+1) * (inow) + snow] = curr_sum + logg[inow];
        istatenow++;
      }
    }
  }
  // convert vector to  3d-array
  IntegerVector dim = {N+1, N+1, days};
  logpolicy.attr("dim") = dim;
  return logpolicy;
}

// [[Rcpp::export]]
List create_fbar_matrix(const double lambda, 
                        const double gamma,
                        const int N,
                        const double c){
  // for every eligible state (st,it), compute the vector storing transition kernel
  // log f(s(t+1), i(t+1) | st,it)
  int ncurrstates = (N+1) * (N+2) / 2;
  int nnextstates = c * c * N + 4 + 4 * c * sqrt(N) + 1 + 1; // this is a crude upper bound
  nnextstates = min(nnextstates , ncurrstates); 
  // Rprintf("matrix dimension is %i by %i\n", ncurrstates, nnextstates);
  // row index correspond to current state
  // column index correspond to next state
  NumericMatrix logfbar(ncurrstates, nnextstates);
  IntegerMatrix nexts(ncurrstates, nnextstates);
  IntegerMatrix nexti(ncurrstates, nnextstates);
  std::fill(nexts.begin(), nexts.end(), -1);
  std::fill(nexti.begin(), nexti.end(), -1);
  std::fill(logfbar.begin(), logfbar.end(), 1.0);
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  int istatenow;
  int cntstatenext;
  istatenow = 0;
  
  for(int snow = N; snow >= 0; snow--){
    for(int inow = N - snow; inow >= 0; inow--){
      // snow and inow defines a current state
      cntstatenext = 0;
      slower = snow * (1 - inow * lambda / N) - c * sqrt(N) / 2;
      supper = snow * (1 - inow * lambda / N) + c * sqrt(N) / 2 + 0.5;
      slower = max(0, slower);
      supper = min(snow, supper);
      for(int snext = supper; snext >= slower; snext--){
        ilower = inow + snow - snext - inow * gamma - c * sqrt(N) / 2;
        iupper = inow + snow - snext - inow * gamma + c * sqrt(N) / 2 + 0.5;
        ilower = max(snow - snext, ilower);
        iupper = min(inow + snow - snext, iupper);
        for(int inext = iupper; inext >= ilower; inext--){
          // snext and inext defines a future state
          logfbar(istatenow, cntstatenext) = sir_logfbar(snow, inow, snext, inext, lambda, gamma, N);
          nexts(istatenow, cntstatenext) = snext;
          nexti(istatenow, cntstatenext) = inext;
          cntstatenext++;
        }
      }
      istatenow++;
      // Rprintf("number of states = %i and analytic nnextstates = %i\n", cntstatenext, nnextstates);
    }
  }
  return List::create(Named("nexts") = nexts,
                      Named("nexti") = nexti,
                      // Named("nextstates") = next_states,
                      Named("logfbar") = logfbar);
}



// assume that logfbar, nexts, nexti already exists
// update these matrices given new values of lambda, gamma
// this function assumes that c is the same as the factor used in create_fbar_matrix
// [[Rcpp::export]]
void update_fbar_matrix(IntegerMatrix & nexts,
                        IntegerMatrix & nexti, 
                        NumericMatrix & logfbar, 
                        const double & lambda,
                        const double & gamma, 
                        const int & N,
                        const double & c){
  // row index correspond to current state
  // column index correspond to next state
  std::fill(nexts.begin(), nexts.end(), -1);
  std::fill(nexti.begin(), nexti.end(), -1);
  std::fill(logfbar.begin(), logfbar.end(), 1.0);
  
  int snow, inow, snext, inext;
  int supper, slower, iupper, ilower;
  int istatenow;
  int cntstatenext;
  istatenow = 0;
  
  for(int snow = N; snow >= 0; snow--){
    for(int inow = N - snow; inow >= 0; inow--){
      // snow and inow defines a current state
      cntstatenext = 0;
      slower = snow * (1 - inow * lambda / N) - c * sqrt(N) / 2;
      supper = snow * (1 - inow * lambda / N) + c * sqrt(N) / 2 + 0.5;
      slower = max(0, slower);
      supper = min(snow, supper);
      for(int snext = supper; snext >= slower; snext--){
        ilower = inow + snow - snext - inow * gamma - c * sqrt(N) / 2;
        iupper = inow + snow - snext - inow * gamma + c * sqrt(N) / 2 + 0.5;
        ilower = max(snow - snext, ilower);
        iupper = min(inow + snow - snext, iupper);
        for(int inext = iupper; inext >= ilower; inext--){
          // snext and inext defines a future state
          logfbar(istatenow, cntstatenext) = sir_logfbar(snow, inow, snext, inext, lambda, gamma, N);
          nexts(istatenow, cntstatenext) = snext;
          nexti(istatenow, cntstatenext) = inext;
          cntstatenext++;
        }
      }
      istatenow++;
    }
  }
}
