#' @title Gibbs sampler for SIS model
#' @description update the agent states for the SIS model given observations such that the complete likelihood is not zero. 
#' It performs N swaps at each time step.
#' @param xx a binary matrix of size model_config$N by length(y), can be initialized from sis_init_hidden_popobservation.R 
#' @param model_config a list containing parameters and network structure
#' @return a new sample
#' 
sis_xx_gibbs_swap_alltimes <- function(xx, model_config){
  if(model_config$network_type == "full"){
    ## we can use the cpp implementation when network is full
    sis_xx_gibbs_swap_full_cpp(xx, model_config$alpha0, model_config$lambda, model_config$gamma);
    return(xx);
  }
  ## choose indices for swap
  ## NOTE: this implementation is O(N)
  ## NOTE: swap moves do not change sum(xt) so the density does not depend on y[t]
  for(time_step in 0 : (dim(xx)[2] - 1)){ 
    ## previous step to now 
    if (time_step == 0){
      alpha_prev2now <- model_config$alpha0;
    }else{
      alpha_prev2now <- sis_get_alpha(xx[,time_step - 1 + 1], model_config);
    }
    ## alpha_prev2now does not change with swaps
    xx_current <- xx[,time_step + 1]
    for(iswap in 1 : dim(xx)[1]){## N swaps at each time
      if (all(xx_current == 0) | all(xx_current == 1)){
        break
      } ## no swap necessary if all 0s or 1s
      i0 <- sample(which(xx_current == 0), size = 1);
      i1 <- sample(which(xx_current == 1), size = 1);
      log_accept_prob <- log(alpha_prev2now[i0]) + log(1 - alpha_prev2now[i1]) - log(alpha_prev2now[i1]) - log(1 - alpha_prev2now[i0]);
      ## now to next step
      if ((time_step + 1) != dim(xx)[2]){
        xx_swap <- xx_current 
        xx_swap[i1] <- 0 
        xx_swap[i0] <- 1
        alpha_now2next_noswap <- sis_get_alpha(xx_current, model_config);
        alpha_now2next_swap <- sis_get_alpha(xx_swap, model_config);
        log_accept_prob <- log_accept_prob + logdbern_sum_cpp(alpha_now2next_swap, xx[, time_step + 1]) - logdbern_sum_cpp(alpha_now2next_noswap, xx[, time_step + 1]);
      }
      ## accept step
      if (log(runif(1)) < log_accept_prob){
        ## swap 
        xx_current[i0] <- 1;
        xx_current[i1] <- 0;
      }
    }
    xx[, time_step + 1] <- xx_current;
  }
  return(xx)
}
