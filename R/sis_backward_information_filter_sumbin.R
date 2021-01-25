#' @title Backward Information Filter for SIS model
#' @description Approximate backward recursion based on coarse-graining the state space.
#' @param y a sequence of observations
#' @param model_config a list containing parameters and network structure
#' @return approximate policy
#' @export 
#' 
sis_backward_information_filter_sumbin <- function(y, model_config){
  ## can implement this function in cpp later
  TT <- length(y)
  full_support_y <- 0 : model_config$N
  ## parameters
  gamma_bar <- mean(model_config$gamma)
  lambda_bar <- mean(model_config$lambda)
  ## approximate transition kernel on i(t+1) | i(t) in a matrix 
  logf_hat <- matrix(-Inf, nrow = 1 +  model_config$N, ncol = 1 + model_config$N)
  logf_hat[1,1] <- 0
  for (it in 1:model_config$N){
    logf_hat[it + 1, ] <- logdsumbin(it, lambda_bar, gamma_bar, model_config$N)
  }
  ## initialize backward information filter from terminal time 
  bif <- matrix(NA, nrow = 1 + model_config$N, ncol = TT)
  bif[,TT] <- dbinom(x = y[TT], size = full_support_y , prob = model_config$rho, log = TRUE)
  ## start the recursion for bif 
  for (tt in (TT - 1): 1){
    ## conditional expectations
    log_cond_expectations <- apply(logf_hat, 1, FUN = function(row) lw.logsum(row + bif[,tt + 1])) 
    ## the line above is the same as log(exp(logf_hat) %*% exp(bif[,tt+1]))
    bif[ , tt] <- dbinom(x = y[tt], size = full_support_y , prob = model_config$rho, log = TRUE) + log_cond_expectations
  }
  return(bif)
}



##### DEBUG CODE
# sis_backward_information_filter_poisbin <- function(y, model_config){
#   ## can implement this function in cpp later
#   TT <- length(y)
#   full_support_y <- 0 : model_config$N
#   ## transition kernel on i(t+1) | i(t) in a matrix 
#   logf <- matrix(-Inf, nrow = 1 +  model_config$N, ncol = 1 + model_config$N)
#   logf[1,1] <- 0
#   for (it in 1:model_config$N){
#     xx <- rep(0,model_config$N)
#     xx[1:it] <- 1
#     a <- sis_get_alpha(agent_state = xx, model_config = model_config)
#     logf[it + 1, ] <- logdpoisbinom_cpp(alpha = a)
#   }
#   ## initialize backward information filter from terminal time 
#   bif <- matrix(NA, nrow = 1 + model_config$N, ncol = TT)
#   bif[,TT] <- dbinom(x = y[TT], size = full_support_y , prob = model_config$rho, log = TRUE)
#   ## start the recursion for bif 
#   for (tt in (TT - 1): 1){
#     ## conditional expectations
#     log_cond_expectations <- apply(logf_hat, 1, FUN = function(row) lw.logsum(row + bif[,tt + 1])) 
#     ## the line above is the same as log(exp(logf_hat) %*% exp(bif[,tt+1]))
#     bif[ , tt] <- dbinom(x = y[tt], size = full_support_y , prob = model_config$rho, log = TRUE) + log_cond_expectations
#   }
#   return(bif)
# }
