#' @title SIR process
#' @description Simulate SIR process as described in Section 4.1
#' @param days length of simulation
#' @param model_config a list containing model parameters, passing the check_model_config function.
#' @param dt increment in time
#' @return a list containing
#' \itemize{
#' \item agent_states : a binary matrix containing the hidden states, size is N by (days + 1)
#' \item y : a vector of reports, length is days + 1
#' } 
#' @export
#' 

sir_simulate <- function(days, model_config, dt = 1){
  agent_states <- matrix(0, nrow = model_config$N, ncol = 1 + days)
  ## simulate initial state according to beta_0
  agent_states[,1] <- (runif(N) < model_config$alpha0)
  if(dt ==1){
    for (d in c(2:(days+1))){
      agent_states[, d] <- sir_kernel(agent_states[, d - 1], model_config); 

    }
  }else{
    for (d in c(2: (days + 1))){
      agent_states[,d] <- sir_simulate_continuous(agent_states[, d - 1], model_config, dt)
    }
  }
  ## partial observation of column sums
  y <- apply(agent_states, 2, function(xx) rbinom(n = 1, size = sum(xx == 1), prob = model_config$rho) )
  return(list(y = y, agent_states = agent_states))
}
