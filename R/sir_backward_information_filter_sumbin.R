#' @title Controlled SMC for SIR model 
#' @description computes the backward information filter
#' @param y population observations
#' @param beta infection rate
#' @param gamma recovery rate
#' @param rho reporting rate
#' @return policy approximate dynamic programming matrix pre-calculated to run an smc filter. This is an array of size [N + 1, N + 1, T + 1]
#' @export
#' 
sir_backward_information_filter_sumbin <- function(y, model_config){
  ## policy[,,t] approximates p(y[t:T]|xt) on log scale
  ## first create the approximate Markov transition kernel
  N <- model_config$N
  num_observations <- length(y)
  policy <- array(data = -Inf, dim = c(model_config$N + 1, model_config$N + 1, num_observations))
  average_lambda <- mean(model_config$lambda)
  average_gamma <- mean(model_config$gamma)
  ## rewrite the commented code in rcpp to save time from the 4 layrs of for-loops
  log_fhat <- function(infection_prev, recovery_prev, infection,recovery){
    if (infection_prev + recovery_prev > N) return(-Inf);
    if (infection + recovery > N ) return(-Inf);
    if (recovery < recovery_prev) return(-Inf);
    if (recovery - recovery_prev > infection_prev) return(-Inf);
    if (N - infection_prev - recovery_prev <  infection - infection_prev + recovery - recovery_prev) return(-Inf);
    dbinom(x = recovery - recovery_prev, size = infection_prev, prob = average_gamma, log = T) + dbinom(x = infection - infection_prev + recovery - recovery_prev, size = N - infection_prev - recovery_prev, prob = average_lambda * infection_prev / N, log = T)
  }
  log_fhat_array <- array(data = -Inf, dim = rep(N + 1,4))
  for (infection_prev in 0:N){
    for (recovery_prev in 0: (N - infection_prev)){
      for (infection in 0:N){
        for (recovery in recovery_prev : N){
          log_fhat_array[infection_prev + 1, recovery_prev + 1, infection + 1, recovery + 1] <- log_fhat(infection_prev , recovery_prev, infection, recovery);
        }
      }
    }
  }
  ## the approximate dynamic programming runs backwards 
  current_logpolicy <- matrix(- Inf, nrow = N + 1, ncol = N + 1)
  for (infection in y[num_observations] : N){
    for (recovery in 0 : (N - infection)){
      current_logpolicy[infection + 1, recovery + 1] <- dbinom(x = y[num_observations], size = infection, prob = model_config$rho, log = T);
    }
  }
  policy[, , num_observations] <- current_logpolicy ;
  for (t in (num_observations - 1) : 1){## dp starts, backwards
    support <- y[t]:N;
    supportsize <- length(support);
    loggt <- dbinom(x = y[t], size = support, prob = model_config$rho, log = T);
    current_logpolicy <- matrix(-Inf, nrow = N + 1, ncol = N + 1);
    for (i in 1 : supportsize){
      infection_prev <- support[i];
      for (recovery_prev in 0 : (N - infection_prev)){
        current_logpolicy[infection_prev + 1, recovery_prev + 1] <- loggt[i] + lw.logsum(policy[ , , t + 1] + log_fhat_array[infection_prev + 1, recovery_prev +1, , ]);
      }
    }
    policy[, , t] <- current_logpolicy;
  }
  return(policy)
}
