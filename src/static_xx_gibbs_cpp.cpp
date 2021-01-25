#include <Rcpp.h>
#include <R.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector static_xx_gibbs_cpp(LogicalVector xx_previous,
                                  NumericVector alpha,
                                  double rho,
                                  int y) {
  int N = xx_previous.length();
  LogicalVector xx(N);
  int sum_xx_previous = sum(xx_previous) ;
  int sum_xx = 0;
  for (int n = 0; n < N; n++){
    sum_xx_previous -= xx_previous[n];
    if (sum_xx_previous + sum_xx < y){
      xx[n] = 1;
    }else{
      double p1 = alpha[n] * R::dbinom(y, 1 + sum_xx + sum_xx_previous, rho, false);
      double p0 = (1 - alpha[n]) * R::dbinom(y, sum_xx + sum_xx_previous, rho, false);
      // sample a uniform random variable 
      GetRNGstate();
      double u = runif(1)(0);
      PutRNGstate();
      xx[n] =  (u < p1 / (p1 + p0));
      // if (n == 0) Rprintf("acceptance probability is %.4f \n", p1 / (p1 + p0)); // debug print message
    }
    sum_xx += xx[n];
  }
  return xx;
}
