#' @title SIS process
#' @description  Markov kernel for SIS process when the increment in time dt < 1
#' @param agent_state binary vector
#' @param model_config
#' @param dt increment in time , dt < 1
sis_simulate_continuous <- function(agent_state, model_config, dt = 0.5){
  steps <- floor(1/dt)
  # dt <- 1/steps
  for (d in 1:steps){
    p <- sis_get_alpha(agent_state, model_config)
    agent_state <- (runif(length(agent_state)) < p)
  }
  return(agent_state)
}
