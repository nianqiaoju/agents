#' @title SIR process
#' @description  Markov kernel for SIR process when the increment in time dt < 1
#' @param agent_state $xt$
#' @param model_config a list containing model parameters, must pass the test of check_model_config.
#' @param dt increment in time , dt < 1
#' @return agent states.
#' 
sir_simulate_continuous <- function(agent_state, model_config, dt = 0.5){
  steps <- floor(1/dt)
  # dt <- 1/steps
  for (d in 1:steps){
    agent_state <- sir_kernel(agent_state, model_config);
  }
  return(agent_state)
}
