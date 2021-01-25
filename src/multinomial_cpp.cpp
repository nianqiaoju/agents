#include <Rcpp.h>
#include "logsumexp.h"
#include "sis_cpp.h"
using namespace Rcpp;
using namespace std;


// zero-based multinomial sampling

// [[Rcpp::export]]
int multinomial_cpp(const NumericVector &logweights, double uniform){
  int n = logweights.size();
  // exponentiate log weights  // not normalized
  double maxlogweights = Rcpp::max(logweights);
  NumericVector cumsum_w = exp(logweights - maxlogweights);
  int draw = 0;
  for (int i = 1; i < n; i++){
    cumsum_w(i) = cumsum_w(i) + cumsum_w(i-1);
  }
  uniform = uniform * cumsum_w(n-1);
  int running_index = 0;
  double sumw = cumsum_w(0);
  if (uniform <= sumw){
    draw = running_index;
  } else {
    while (uniform > sumw){
      sumw = cumsum_w(running_index);
      running_index ++;
    }
    draw = running_index-1;
  }
  return draw;
}