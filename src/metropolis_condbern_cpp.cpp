#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// This function gives one approximate sample from the Conditional Bernoulli distribution
// specified by an integer sum_x and a vector of probabilities alpha
// The MCMC chain is initialized with sampling uniform without replacements
// it uses the swap kernel
// it runs for num_mcmc iterations

// [[Rcpp::export]]
LogicalVector metropolis_condbern_cpp(int sum_x,
                                      NumericVector alpha,
                                      int num_mcmc) {
  int N_ = alpha.length();
  LogicalVector xx(N_); //
  IntegerVector n1(sum_x); // locations of 1s
  IntegerVector n0(N_ - sum_x); // locations of 0s
  NumericVector logw = log(alpha) - log(1 - alpha); // log odds
  // Rcout << "odds = " << logw << "\n";

  // instantiate the chain 
  n1 = sample(N_, sum_x, false); // sample without replacement
  IntegerVector n01 = seq_len(N_);
  n0 = setdiff(n01, n1);
  // Rcout << "n1 = " << n1 << "\n";
  // Rcout << "n0 = " << n0 << "\n";

  // swap moves
  int j0; // index in n0, starts at 0
  int j1; // index in n1, starts at 1
  double log_accept; // log acceptance probability

  for(int iter = 0; iter < num_mcmc; iter++){
  	  GetRNGstate();
      NumericVector uniforms = runif(3); 
      PutRNGstate();
      j1 = floor(uniforms[0] * sum_x);
      j0 = floor(uniforms[1] * (N_ - sum_x));
      log_accept = logw[n0[j0] - 1] - logw[n1[j1] - 1];
      if(log(uniforms[2]) < log_accept){
      	int temp = n1[j1];
      	n1[j1] = n0[j0];
      	n0[j0] = temp;
      }
  }
  for(int j1 = 0; j1 < sum_x; j1++){
  	xx[n1[j1] - 1] = 1;
  }
  return(xx);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
sum_x <- 3
alpha <- runif(10)
num_mcmc <- 20
metropolis_condbern_cpp(sum_x, alpha, num_mcmc)
*/
