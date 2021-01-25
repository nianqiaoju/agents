#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// This function gives one exact sample from the Conditional Bernoulli distribution
// using the idchecking method from Chen and Liu 1997

// [[Rcpp::export]]
LogicalVector idchecking_cpp(const int sum_x,
                   const NumericVector &alpha,
                   const NumericVector &random_number) {
  int N = alpha.length();
  NumericMatrix logq(N+1,N);
  std::fill(logq.begin(), logq.end(), R_NegInf);
  // compute all the intermediate q(h,n) in log scale 
  // q(h,n) = PB(h; alpha[n : N])
  // base cases: 
  // (1) q(0,n) = prod(1 - alpha[n:N]) for all n 
  // (2) q(h, N - h + 1) = prod(alpha[N - h + 1 : N])
  // (3) q(h, n ) = 0 for h > N - n + 1
  // recursion: q(h,n) = alpha[n] * q(h-1, n+1) + (1 - alpha[n]) * q(h, n+1)
  // logq[h , n] = log(q(h, n)) 
  logq(0,N - 1) = log(1 - alpha[N - 1]);
  for(int n = N - 1; n > 0; n--){
    logq(0,n - 1) = logq(0, n - 1 + 1) + log(1 - alpha[n - 1]);
  }
  logq(1,N - 1) = log(alpha[N - 1]);
  double ls1, ls2, maxls;
  for(int n = N - 1; n >= 1; n--){
    for(int h = 1; (h <= sum_x) && (h <= N - n + 1); h++){
      ls1 = log(alpha[n - 1]) + logq(h - 1, n + 1 - 1);
      ls2 = log(1 - alpha[n - 1]) + logq(h , n + 1 - 1);
      maxls = max(ls1, ls2);
      logq(h , n - 1) = log(exp(ls1 - maxls) + exp(ls2 - maxls)) + maxls;
      if(Rcpp::traits::is_nan<REALSXP>(logq(h, n-1))){
        logq(h, n - 1) = R_NegInf;
      }
    }
  }
  // Rcout << "The value of logq : " << logq << "\n";

  LogicalVector xx(N);
  int h = 0; // partial sum
  double logp = R_NegInf;
  // probablity of xx[n] = 1 is alpha[n] * q(i - h - 1, n + 1) / q(i - h, n)
  for(int n = 1; n < N ; n++){
    logp = log(alpha[n - 1]) + logq(sum_x - h - 1, n + 1 - 1) - logq(sum_x - h, n - 1);
    if(log(random_number[n - 1]) < logp){
      xx[n - 1] = true;
      h++;
    }
    else{
      xx[n - 1] = false;
    }
    if(h == sum_x) break;
  }
  xx[N - 1] = (h == sum_x - 1);
  return(xx);
}