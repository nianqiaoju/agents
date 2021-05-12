#' @title transition probabilty of the SIS process
#' @description for the sis model, given xt and static attributes, computes the probability of x[t+1,k] = 1 for each k 
#' @param agent_state the state at time t 
#' @param model_config a list of model parameters and network structure
#' @return the alpha vector
#' @export  

sis_get_alpha <- function(agent_state, model_config){
  ## assumes that agent_state is a vector
  if(model_config$network_type == 'full'){
    alpha <- sis_get_alpha_full_cpp(agent_state, model_config$lambda, model_config$gamma);
  }else{
    alpha <- agent_state * ( 1 - model_config$gamma) + (1 - agent_state) * (model_config$lambda * ( agent_state %*% model_config$adjacency_matrix_b));
  }
  return(as.vector(alpha))
}