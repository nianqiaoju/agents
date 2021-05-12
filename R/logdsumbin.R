#' @title Sum of two independent Binomials 
#' @description computes the log density of SumBinomial distribition, parametrized by lbar, gbar and it.
#' @param it number of infections.
#' @param lbar average infection rate.
#' @param gbar average recovery rate.
#' @param N integer, must satisfy N = length(lbar).
#' @return log transition densities p(i(t+1) = 0 : N | it, lbar, gbar);

logdsumbin <- function(it,lbar,gbar, N){
  new_infection <- dbinom(x = 0 : (N - it), size = N - it, prob = lbar * it / N, log = TRUE);
  stay_infection <- dbinom(x = 0: it, size = it, prob = 1 - gbar, log = TRUE);
  logdtransitions <- rep(-Inf, 1 + N);
  ## all vectors are offset by 1 since index start at zero.
  ## convolution: logdtransitions[z + 1] = sum of (new_infection[y + 1] + stay_infection[z - y + 1]) for y = 0 : z.
  for(new_i in 0: (N - it)){
    for(stay_i in 0 : it){
      a <- logdtransitions[new_i + stay_i + 1];
      b <- new_infection[new_i + 1] + stay_infection[stay_i + 1];
      mab <- max(a,b);
      a <- a - mab;
      b <- b - mab;
      logdtransitions[new_i + stay_i + 1] <- log(exp(a) + exp(b)) + mab;
    }
  }
  return(logdtransitions);
}