#' @title Controlled SMC for SIR model 
#' @description computes the backward information filter as a function of s and i. This function will replace the current sir_backward_information_filter_sumbin once completed.
#' @param y population observations
#' @param  model_config a list containing model parameters and network structure
#' @return policy approximate dynamic programming matrix pre-calculated to run an smc filter. This is an array of size [N + 1, N + 1, T + 1]
#' @export
#' 
sir_backward_information_filter_si <- function(y, model_config){
  ## policy[,,t] approximates p(y[t:T]|xt) on log scale
  ## first create the approximate Markov transition kernel
  N <- model_config$N;
  num_observations <- length(y);
  logpolicy <- array(data = -Inf, dim = c(model_config$N + 1, model_config$N + 1, num_observations));
  ## logpolicy[s+1,i+1,t+1] is log(psi_t(s,i)).
  average_lambda <- mean(model_config$lambda);
  average_gamma <- mean(model_config$gamma);
  # rewrite line 18 - 34 in rcpp to save time from the 4 layrs of for-loops
  log_fhat <- function(scount_prev, icount_prev, scount, icount){
    if (icount_prev + scount_prev > N) return(-Inf);
    if (icount + scount > N ) return(-Inf);
    if (scount > scount_prev) return(-Inf);
    if (icount_prev + scount_prev < icount + scount ) return(-Inf);
    dbinom(x = scount_prev + icount_prev - icount - scount, size = icount_prev, prob = average_gamma, log = T) + dbinom(x = - scount + scount_prev, size = scount_prev, prob = average_lambda * icount_prev / N, log = T)
  }
  log_fhat_array <- array(data = -Inf, dim = rep(N + 1,4))
  for (icount_prev in 0:N){
    for (scount_prev in 0: (N - icount_prev)){
      for (icount in 0:N){
        for (scount in 0 : scount_prev){
          log_fhat_array[scount_prev + 1, icount_prev + 1, scount + 1,  icount + 1] <- log_fhat(scount_prev , icount_prev, scount, icount);
        }
      }
    }
  }
  ## the approximate dynamic programming runs backwards 
  ## terminal time 
  current_logpolicy <- matrix(- Inf, nrow = N + 1, ncol = N + 1)
  for (icount in y[num_observations] : N){
    for (scount in 0 : (N - icount)){
      current_logpolicy[scount + 1, icount + 1] <- dbinom(x = y[num_observations], size = icount, prob = model_config$rho, log = T);
    }
  }
  logpolicy[, , num_observations] <- current_logpolicy ;
  ## t <= terminal time - 1;
  for (t in (num_observations - 1) : 1){## dynamic programming starts, backwards
    support <- y[t]:N;
    supportsize <- length(support);
    loggt <- dbinom(x = y[t], size = support, prob = model_config$rho, log = T);
    current_logpolicy <- matrix(-Inf, nrow = N + 1, ncol = N + 1);
    for (i in 1 : supportsize){
      icount_prev <- support[i];
      for (scount_prev in 0 : (N - icount_prev)){
        current_logpolicy[scount_prev + 1, icount_prev + 1] <- loggt[i] + lw.logsum(logpolicy[ , , t + 1] + log_fhat_array[scount_prev + 1, icount_prev + 1, , ]);
      }
    }
    logpolicy[, , t] <- current_logpolicy;
  }
  return(logpolicy)
}
