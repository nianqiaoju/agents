#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix sis_forward_algorithm_cpp(NumericVector logf0, 
                                        NumericMatrix logdtransition,  
                                        NumericVector y, 
                                        NumericVector all_sum_x, 
                                        double rho) {
  int size_state_space = logdtransition.nrow(); // this number is in fact 2**N for SIS or 3**N for SIR
  int num_observations = y.size();
  NumericMatrix m(size_state_space, num_observations); // each row is an individual and each column is a day 
  for (int istate = 0; istate < size_state_space; istate ++){
    m(istate,0) = logf0(istate) + R::dbinom(y(0), all_sum_x(istate), rho, true); 
  }
  for (int t = 1; t <= num_observations - 1; t++){
  	for (int istate = 0; istate <  size_state_space ; istate++){
  		NumericVector logweights = m(_, t - 1) + logdtransition( _ , istate);
  		if(all(is_infinite(logweights))){
  		  m(istate, t) = R_NegInf;
  		}else{
  		  double maxlogweights = max(logweights);
  		  m(istate,t) = maxlogweights + log(sum(exp(logweights - maxlogweights))) + R::dbinom(y(t), all_sum_x(istate), rho, true);
  		}
  	}
  }
  return m;
}



