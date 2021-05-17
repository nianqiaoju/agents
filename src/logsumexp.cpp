#include "logsumexp.h"
using namespace Rcpp;
using namespace std;
//' logsum using cpp 
//' @param lw log weights.
//' @export 
// [[Rcpp::export]]
double lw_logsum_cpp(const NumericVector &lw){
  // input: lw - log weights 
  // output: sum of weights, i.e. log(sum(exp(lw)))
  double maxlw = max(lw);
  if (traits::is_infinite<REALSXP>(maxlw)){
    return R_NegInf;
  }
  else{
    return maxlw + log(sum(exp(lw - maxlw)));
  }
}


//' normalize logweights
//' @param lw log weights.
//' @export 
// [[Rcpp::export]]
NumericVector lw_normalize_cpp(const NumericVector & lw){
  // input: lw - log weights 
  // output: normalized weights, aka they must sum to 1 
  NumericVector weights (lw.size());
  double maxlw = max(lw);
  if (traits::is_infinite<REALSXP>(maxlw)){
    return weights; // return zeros if sum of the weights is 0 
  }
  else{
    NumericVector logweights (lw.size());
    double sum_weights = 0;
    for(int i = 0; i < lw.size(); i++){
      logweights[i] = lw[i] -  maxlw;
      sum_weights += exp(logweights[i]);
    }
    // sum_w is cumsum now
    for(int i = 0; i < lw.size(); i++){
      weights[i] = exp(logweights[i]) / sum_weights;
    }
    return weights;
  }
}

/*
normalize logweights to weights (in place) 
and return the log normalizing constant.
This function is used in sir_apf.
@param lw log weights
@export
 */
// [[Rcpp::export]]
double lw_logsum_normalize_cpp(NumericVector &lw){
  double sum_weights = 0;
  double maxlw = max(lw);
  if(traits::is_infinite<REALSXP>(maxlw)){
    return R_NegInf;
  }else{
    for(int i = 0; i < lw.size(); i++){
      lw[i] = lw[i] - maxlw;
      sum_weights += exp(lw[i]);
    }
    for(int i = 0; i < lw.size(); i++){
      lw[i] = exp(lw[i]) / sum_weights;
    }
    return maxlw + log(sum(exp(lw - maxlw)));
  }
}