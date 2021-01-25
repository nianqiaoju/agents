#' @title Initalize Gibbs sampler for SIS model
#' @description initialize the agent states for the SIS model given observations such that the complete likelihood is not zero.
#' @param y a list of population observations 
#' @param model_config a list that must contain
#' \itemize{
#' \item N : population size
#' \item rho : reporting rate
#' }
#' @return agent_states

sis_xx_initialize <- function(y, model_config){
  num_observations <- length(y)
  ## instantiate the agent states
  agent_states <- matrix(FALSE, nrow = model_config$N, ncol = num_observations)
  ## sample total number of infections for each day such that 
  ## (1)  yt <= it <= N and (2) 1 <= it
  sum_agent_states <- sapply(y, function(x) x + rbinom(n = 1, size = x, prob = model_config$rho))
  sum_agent_states[is.na(sum_agent_states)] <- 1 ## the case y[t] == 0 
  sum_agent_states <- pmin(sum_agent_states, model_config$N)  
  sum_agent_states <- pmax(sum_agent_states, 1)
  ## randomly assign infections for each day
  for (t in 1:num_observations){
    agent_states[sample.int(n = model_config$N, size = sum_agent_states[t], replace = FALSE) , t] <- TRUE
  } 
  return(agent_states)
}
