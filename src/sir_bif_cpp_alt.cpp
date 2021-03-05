#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include<algorithm>

//

// [[Rcpp::export]]
NumericVector sir_bif_cpp_alt(const IntegerVector y, 
                          const double lambda, 
                          const double gamma, 
                          const double rho, 
                          const int N, 
                          const double c) {
  int days = y.size();
  // Rcout << "observation length = " <<  days <<"\n";
  // initialize a vector to -Inf
  NumericVector logpolicy((N+1) * (N+1)  * days);
  std::fill(logpolicy.begin(), logpolicy.end(), R_NegInf);
  // deal with the terminal time first
  int t = days - 1;
  // Rprintf("at day t = %i\n", t);
  for(int icount = N; icount >= y[t]; icount--){
    for(int scount = 0; scount + icount <= N; scount++){
      // Rprintf("s = %i, i = %i\n", scount, icount);
      logpolicy[t * (N+1) * (N+1) + (N+1) * (icount) + scount] = R::dbinom(y[t], icount, rho, true);
    }
  }
  int slower, supper, ilower, iupper;
  int icount, scount, scount_next, icount_next;
  double logg;
  double curr_sum, add_item, curr_max;
  for(int t = days - 2; t >= 0; t--){
    // Rprintf("at day t = %i\n", t);
    // compute the conditional expecation for each value of st and it
    for(int icount = y[t]; icount <= N; icount++){
      logg = R::dbinom(y[t], icount, rho, true);
      for(int scount = 0; icount + scount <= N; scount++){
        // Rprintf("s = %i, i = %i\n", scount, icount);
        curr_sum = R_NegInf; // initialize the conditonal expectation to 0
        slower = scount * (1 - icount * lambda / N) - c * sqrt(N);
        supper = scount * (1 - icount * lambda / N) + c * sqrt(N);
        slower = max(0, slower);
        supper = min(scount, supper);
        supper = min(N - y[t+1], supper);
        // Rprintf("snext from %i to %i\n", slower, supper);
        // conditonal expectation is summation over scount_next and icount_next
        for(int scount_next = slower; scount_next <= supper; scount_next++){
          ilower = icount + scount - scount_next - icount * gamma - c * sqrt(N);
          iupper = icount + scount - scount_next - icount * gamma + c * sqrt(N);
          ilower = max(icount + scount - scount_next - icount, ilower);
          ilower = max(y[t+1], ilower);
          iupper = min(icount + scount - scount_next, iupper);
          // Rprintf("inext from %i to %i\n", ilower, iupper);
          for(int icount_next = ilower; icount_next <= iupper; icount_next++){
            add_item = R::dbinom(scount + icount - icount_next - scount_next, icount, gamma, true) +
              R::dbinom(scount - scount_next, scount, lambda * icount / N, true) +
              logpolicy[(t + 1) * (N+1) * (N+1) + (N+1) * (icount_next) + scount_next];
            if(!Rcpp::traits::is_infinite<REALSXP>(add_item)){
              curr_max = max(curr_sum, add_item);
              curr_sum = log(exp(curr_sum - curr_max) + exp(add_item - curr_max)) + curr_max;
            }
          }
        }
        // psi_t(st,it) is observation density * conditional expectation
        logpolicy[t * (N+1) * (N+1) + (N+1) * (icount) + scount] = curr_sum + logg;
      }
    }
  }
  // convert vector to  3d-array
  IntegerVector dim = {N+1, N+1, days};
  logpolicy.attr("dim") = dim;
  return logpolicy;
}

