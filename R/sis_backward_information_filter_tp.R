#' @title Backward information filter for SIS model
#' @description This function computes performs the backward information filters steps. 
#' It returns a matrix of size (N + 1) by T and each column approximates observation densities of 
#' observing current and all the future observations given the current state i(xt).
#' @param y population observations
#' @param model_config a list that contains parameters and network structure
#' @return policy approximate dynamic programming matrix pre-calculated to run an smc filter
#' @export

sis_backward_information_filter_tp <- function(y, model_config){
  num_observations <- length(y)
  policy <- matrix(NA, nrow = model_config$N + 1, ncol = num_observations) ## index + 1 for starting at 0
  average_lambda <- mean(model_config$lambda)
  average_gamma <- mean(model_config$gamma)
  ## initialize the backward recursion at t = T
  current_logpolicy <- dbinom(y[num_observations], 
                              size = 0 : model_config$N, 
                              prob = model_config$rho, 
                              log = TRUE);
  policy[, num_observations] <- current_logpolicy;
  current_support <- y[num_observations]:model_config$N; ## bif is not -Inf only at these values of it
  for (t in (num_observations-1):1){## dp starts, backwards
    support <- y[t]:N 
    ## for each t
    ## support contains y[t] : N;
    ## current_support contains y[t+1] : N
    support_size <- length(support)
    logcondexp <- rep(-Inf, support_size)
    for (support_index in 1 : support_size){
      i <- support[support_index]
      # compute mean of translated Poisson approximation
      approx_mean <- (model_config$N - i) * (average_lambda * i / N) + i * (1 - average_gamma)
      # compute variance of translated Poisson approximation
      approx_variance <- approx_mean - (N - i) * (average_lambda * i / N)^2 - (1 - average_gamma)^2 * i
      # evaluate translated Poisson PMF
      # these are probabilities of p(i(t+1) | it)
      transpoi_logpmf <- rep(-Inf, N + 1) ## the plus 1 is indexing for starting at zero.
      transpoi_logpmf[current_support + 1] <- logdtranspoisson(current_support, approx_mean, approx_variance)
      # compute approximate conditional expectation
      logsummand <- transpoi_logpmf + current_logpolicy
      if (all(is.infinite(logsummand))){
        logcondexp[support_index] <- -Inf
      } else{
        maxlogsummand <- max(logsummand)
        summand <- exp(logsummand - maxlogsummand)
        logcondexp[support_index] <- log(sum(summand)) + maxlogsummand
      }
    }
    # compute new policy
    current_logpolicy <- dbinom(y[t], 
                                size = 0:model_config$N, 
                                prob = model_config$rho, 
                                log = TRUE) ## this is the binomial part  
    current_logpolicy[support + 1] <- current_logpolicy[support + 1] + logcondexp 
    policy[, t] <- current_logpolicy
    current_support <- support
  }
  return(policy)
}
