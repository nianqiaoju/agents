#include <Rcpp.h>
#include "logdensities.h"
using namespace Rcpp;
using namespace std;


// alpha is a length N vector of success probabilities 
// the poisson binomial distribution takes value from 0 to N
// returns the probability mass function as a length N + 1 vector
  
// [[Rcpp::export]]
NumericVector logdpoisbinom_cpp(const NumericVector &alpha) {
  int N = alpha.length();
  NumericVector previous_log_densities (N + 1, R_NegInf);
  NumericVector current_log_densities (N + 1, R_NegInf);
  
  // initilaize the recursion, start with q(0,N) and q(1,N);
  previous_log_densities[0] = log(1 - alpha[N - 1]);
  previous_log_densities[1] = log(alpha[N - 1]);
  
  double ls1, ls2, maxls;
  for(int n = N - 1; n > 0; n--){
    current_log_densities[0] = previous_log_densities[0] + log(1 - alpha[n - 1]);
    for(int h = 1; h <= N - n + 1 ; h++){
      ls1 = log(alpha[n - 1]) + previous_log_densities[h - 1];
      ls2 = log(1 - alpha[n - 1]) + previous_log_densities[h];
      maxls = max(ls1, ls2);
      current_log_densities[h] = log(exp(ls1 - maxls) + exp(ls2 - maxls)) + maxls;
    }
    // assign current_log_densities to previous_log_densities
    // previous_log_densities = current_log_densities;
    for(int h = 0; h <= N; h++){
      // replace na with -inf
      if (Rcpp::traits::is_nan<REALSXP>(current_log_densities[h])){
        previous_log_densities[h] = R_NegInf;
      }else{
        previous_log_densities[h] = current_log_densities[h];
      }
    }
    std::fill(current_log_densities.begin(), current_log_densities.end(), R_NegInf);
  }
  return previous_log_densities;
}


// x is a binary vector such that each entry is 
// independent Bernoulli with success probability alpha
// returns the log joint probability mass function


// [[Rcpp::export]]
double logdbern_sum_cpp (const NumericVector & alpha, const  LogicalVector &x){
  // if alpha is zero vector and x is zero vector, then probability is one
  double lprod = 0;
  int N = x.size();
  for (int n = 0; n < N; n++){
    lprod += ((x[n] == 1)? log(alpha[n]) : log(1 - alpha[n]));
  }
  return(lprod);
}

/*** R
alpha <- c(1,0,0.5)
logdpoisbinom_cpp(alpha)
logdpoisbinom(alpha)
*/
