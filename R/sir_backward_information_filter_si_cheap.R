#' @title Controlled SMC for SIR model 
#' @description computes the backward information filter as a function of s and i. This function will replace the current sir_backward_information_filter_sumbin once completed.
#' @param y population observations
#' @param  model_config a list containing model parameters and network structure
#' @return policy approximate dynamic programming matrix pre-calculated to run an smc filter. This is an array of size [N + 1, N + 1, T + 1]
#' @export
#' 
sir_backward_information_filter_si_cheap <- function(y, model_config, c = 2){
  ## policy[,,t] approximates p(y[t:T]|xt) on log scale
  ## first create the approximate Markov transition kernel
  N <- model_config$N;
  num_observations <- length(y);
  logpolicy <- array(data = -Inf, dim = c(model_config$N + 1, model_config$N + 1, num_observations));
  ## logpolicy[s+1,i+1,t+1] is log(psi_t(s,i)).
  average_lambda <- mean(model_config$lambda);
  average_gamma <- mean(model_config$gamma);
  logfbar <- function(scount, icount, scount_next, icount_next){
    dbinom(x = scount + icount - icount_next - scount_next, size = icount, prob = average_gamma, log = T) + dbinom(x = scount - scount_next, size = scount, prob = average_lambda * icount / N, log = T);
  }
  # log_fhat_array <- array(data = -Inf, dim = rep(N + 1,4))
  # for (icount in 0:N){
  #   for (scount in 0: (N - icount_prev)){
  #     for (scount_next in 0: scount_prev){
  #       for (scount_next in 0 : scount_prev){
  #         log_fhat_array[scount_prev + 1, icount_prev + 1, scount + 1,  icount + 1] <- log_fhat(scount_prev , icount_prev, scount, icount);
  #       }
  #     }
  #   }
  # }
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
      icount <- support[i];
      for (scount in 0 : (N - icount)){
        curr_sum <- - Inf;
        ## line 50-52 correspond to the original function
        # for(scount_next in 0 : scount){
        #   for(icount_next in (scount - scount_next) : (icount + scount - scount_next)){
        for(scount_next in max(0, floor(scount * (1 - icount * average_lambda / N) - c * sqrt(scount * (1 - icount * average_lambda / N) * icount * average_lambda / N))) :  min(scount, floor(scount * (1 - icount * average_lambda / N) + c * sqrt(scount * (1 - icount * average_lambda / N) * icount * average_lambda / N)))){
          for(icount_next in (icount + scount - scount_next - min(icount , floor(icount * average_gamma + c * sqrt(icount * average_gamma * (1 - average_gamma))))) : (icount + scount - scount_next - floor(max(0, icount * average_gamma - c * sqrt(icount * average_gamma * (1 - average_gamma)))))){
            ## add the value of logfbar(scount_next, icount_next | scount, icount ) + logpolicy(scount, icount, t)  
            add_item <- logfbar(scount, icount, scount_next, icount_next) + logpolicy[scount_next + 1, icount_next + 1,t + 1];
            curr_max <- max(curr_sum, add_item);
            if(is.finite(curr_max)){
              curr_sum <- curr_sum - curr_max;
              add_item <- add_item - curr_max;
              curr_sum <- log(exp(curr_sum) + exp(add_item)) + curr_max;
            }
          }
        }
        current_logpolicy[scount + 1, icount + 1] <- loggt[i] + curr_sum;
      }
    }
    logpolicy[, , t] <- current_logpolicy;
  }
  return(logpolicy)
}
