#' @title Controlled SMC for SIR model 
#' @description computes the backward information filter
#' @param y population observations
#' @param beta infection rate
#' @param gamma recovery rate
#' @param rho reporting rate
#' @return policy approximate dynamic programming matrix pre-calculated to run an smc filter. This is an array of size [N + 1, N + 1, T + 1]
#' @export
#' 
### this function is under-development and will replace sir_backward_information_filter_sumbin once tested.
sir_backward_information_filter_sumbin_dev <- function(y, model_config, index2i, index2r){
  ## policy[,,t] approximates p(y[t:T]|xt) on log scale
  ## first create the approximate Markov transition kernel
  N <- model_config$N;
  num_observations <- length(y);
  irspace_size <- (model_config$N + 1) * (model_config$N + 2) / 2
  policy_array <- array(data = -Inf, dim = c(model_config$N + 1, model_config$N + 1, num_observations))
  policy_matrix <- matrix(NA, nrow = irspace_size, ncol = num_observations);
  average_lambda <- mean(model_config$lambda)
  average_gamma <- mean(model_config$gamma)
  ## compute log policy
  policy_matrix <- sir_bif_policy_matrix_cpp(average_lambda, average_gamma, y, N, model_config$rho);
  ## convert the policy_matrix to policy_array, to be used in cSMC
  for(t in 1 : num_observations){
    for(ir_index in 1 : irspace_size){
      policy_array[index2i[ir_index] + 1, index2r[ir_index] + 1,t] <- policy_matrix[ir_index,t];
    }
  }
  return(policy_array)
}

##  below is the R implementation of sir_bif_policy_cpp
## use vectors to track index2i and index2r
# index2i <- index2r <- integer(irspace_size);
# for (i in 0 : model_config$N){
#   for (r in 0 : (model_config$N - i)){
#     ir2index <- misc_ir2index_cpp(i,r,N);
#     index2i[ir2index + 1] <- i ;
#     index2r[ir2index + 1] <- r;
#   }
# }
# log_fhat <- function(infection_prev, recovery_prev, infection,recovery){
#   if (infection_prev + recovery_prev > N) return(-Inf);
#   if (infection + recovery > N ) return(-Inf);
#   if (recovery < recovery_prev) return(-Inf);
#   if (recovery - recovery_prev > infection_prev) return(-Inf);
#   if (N - infection_prev - recovery_prev <  infection - infection_prev + recovery - recovery_prev) return(-Inf);
#   dbinom(x = recovery - recovery_prev, size = infection_prev, prob = average_gamma, log = T) + dbinom(x = infection - infection_prev + recovery - recovery_prev, size = N - infection_prev - recovery_prev, prob = average_lambda * infection_prev / N, log = T)
# }
# log_fhat_matrix <-  matrix(NA, nrow = irspace_size, ncol = irspace_size);
# for(prev_state in 1 : irspace_size){
#   for(next_state in 1 : irspace_size){
#     log_fhat_matrix[prev_state, next_state] <- log_fhat(index2i[prev_state],index2r[prev_state], index2i[next_state], index2r[next_state]);
#   }
# }
# ## the approximate dynamic programming runs backwards 
# current_logpolicy <- rep(-Inf, irspace_size);
# for (infection in y[num_observations] : N){
#   for (recovery in 0 : (N - infection)){
#     current_logpolicy[misc_ir2index_cpp(infection, recovery, N) + 1] <- dbinom(x = y[num_observations], size = infection, prob = model_config$rho, log = T);
#   }
# }
# policy_matrix[,num_observations] <- current_logpolicy;
# for (t in (num_observations - 1) : 1){## dp starts, backwards
#   support <- y[t]:N;
#   supportsize <- length(support);
#   loggt <- dbinom(x = y[t], size = support, prob = model_config$rho, log = T);
#   current_logpolicy <- rep(-Inf, irspace_size);
#   for (i in 1 : supportsize){
#     infection_prev <- support[i];
#     for (recovery_prev in 0 : (N - infection_prev)){
#       current_logpolicy[misc_ir2index_cpp(infection_prev, recovery_prev, N) + 1]  <- loggt[i] + lw.logsum(policy_matrix[ , t + 1] + log_fhat_matrix[misc_ir2index_cpp(infection_prev,recovery_prev,N) + 1,]);
#     }
#   }
#   policy_matrix[, t] <- current_logpolicy;
# }