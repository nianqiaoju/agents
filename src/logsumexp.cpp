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
for each column, normalize logweights to weights (in place) 
and return the log normalizing constants in a vector.
This function is used in sir_apf.
@param lw a matrix, each column is a vector of log weights
@return a vector containing the log normalizing constant of each column
@export
 */
// [[Rcpp::export]]
NumericVector lw_logsum_normalize_byCol(NumericMatrix &lw){
  double sum_weights = 0;
  double maxlw;
  NumericVector logsums(lw.ncol());
  for(int icol = 0; icol < lw.ncol(); icol++){
    maxlw = max(lw( _, icol));
    sum_weights = 0;
    if(traits::is_infinite<REALSXP>(maxlw)){
      logsums[icol] = R_NegInf;
    }else{
      for(int irow = 0; irow < lw.nrow(); irow++){
        lw(irow,icol) = lw(irow, icol) - maxlw;
        sum_weights += exp(lw(irow,icol));
      }
      for(int irow = 0; irow < lw.nrow(); irow++){
        lw(irow, icol) = exp(lw(irow, icol)) / sum_weights; // lw has become the normalized weights
      }
      logsums[icol] = maxlw + log(sum_weights);
    }
  }
  return logsums;
}