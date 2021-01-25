#' @title Sum of two independent Binomials 
#' @description computes the log density of SumBinomial distribition, parametrized by lambda_bar, gamma_bar and it
#' @param it number of infections
#' @param lambda_bar average infection rate
#' @param gamma_bar average recovery rate
#' @return log transition densities p(i(t+1) = 0 : N | it, lambda_bar, gamma_bar);

logdsumbin <- function(i_current,lambda_bar,gamma_bar, N){
  new_infection <- dbinom(x = 0 : (N - i_current), size = N - i_current, prob = lambda_bar * i_current / N, log = TRUE)
  stay_infection <- dbinom(x = 0: i_current, size = i_current, prob = 1 - gamma_bar, log = TRUE)
  logdtransitions <- rep(-Inf, 1 + N)
  for(new_i in 0: (N - i_current)){
    for(stay_i in 0 : i_current){
      ab <- rep(NA, 2)
      ab[1] <- logdtransitions[new_i + stay_i + 1]
      ab[2] <- new_infection[new_i + 1] + stay_infection[stay_i + 1]
      logdtransitions[new_i + stay_i + 1] <- lw.logsum(ab)
    }
  }
  return(logdtransitions)
}