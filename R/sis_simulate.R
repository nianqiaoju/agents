#' @title SIS process
#' @description Simulate SIS process as described in Section 3.1
#' @param days length of simulation
#' @param model_config a list containing
#' \itemize{
#' \item N
#' \item adjacency matrix
#' \item alpha0
#' \item gamma
#' \item beta
#' \item rho}
#' @param dt increment in time
#' @return a list containing
#' \itemize{
#' \item agent_states : a binary matrix containing the hidden states, size is N by (days + 1)
#' \item y : a vector of reports, length is days + 1
#' \item 
#' } 
#' @export
#' 

sis_simulate <- function(days, model_config, dt = 1){
  agent_states <- matrix(0, nrow = model_config$N, ncol = 1 + days)
  ## simulate initial state according to beta_0
  agent_states[,1] <- (runif(N) < model_config$alpha0)
  if(dt ==1){
    for (d in c(2:(days+1))){
      ##  need the function state_to_alpha
      alpha_t <- rep(NA, N)
      alpha_t <- sis_get_alpha(agent_states[,d - 1], model_config)
      agent_states[,d] <- (runif(N) < alpha_t)
    }
  }else{
    for (d in c(2: (days + 1))){
      agent_states[,d] <- sis_simulate_continuous(agent_states[, d - 1], model_config, dt)
    }
  }
  ## partial observation of column sums
  y <- apply(agent_states, 2, function(xx) rbinom(n = 1, size = sum(xx), prob = model_config$rho) )
  return(list(agent_states = agent_states, y = y))
}


