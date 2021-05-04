#' @title SIR model transition kernels
#' @rdname sir_get_alpha
#' @description the alpha function for SIR process.
#' @param agent_states the agent states at some time t 
#' @param model_config a list containing model parameters
#' @return \code{sir_get_alpha_sir(agent_states,model_config)} returns a matrix with three rows : alpha_s, alpha_i and alpha_r.
#' @export
#' 
sir_get_alpha_sir <- function(agent_states, model_config){
	alpha_matrix <- matrix(NA, nrow = 3, ncol = model_config$N);
	## probability of recovery
	alpha_matrix[3,] <- (agent_states == 2) + (agent_states ==1 ) * model_config$gamma
	alpha_matrix[1,] <- (agent_states == 0) * pmax(1e-4, 1 - model_config$lambda * t(agent_states == 1) %*% model_config$adjacency_matrix_b)
	# alpha_matrix[1,] <- (agent_states == 0) * (1 - model_config$lambda * t(agent_states == 1) %*% model_config$adjacency_matrix_b)
	alpha_matrix[2, ] <- 1 - alpha_matrix[3,] - alpha_matrix[1,]
	return(alpha_matrix)
}

#' @rdname sir_get_alpha
#' @return \code{sir_get_alpha_i(agent_states, model_config)} returns alpha_i
#' @export
sir_get_alpha_i <- function(agent_states, model_config){
  if(model_config$network_type == "full"){
    it <- sum(agent_states == 1);
    alpha <- (agent_states == 0) * pmax(1e-4, model_config$lambda * it / N) + (agent_states == 1) * (1 - model_config$gamma);
  }else{
    alpha <- (agent_states == 0) * pmax(1e-4, model_config$lambda * (agent_states == 1) %*% model_config$adjacency_matrix_b)  + (agent_states == 1) * ( 1 - model_config$gamma);
    
  }
	return(alpha)
}


#' @rdname  sir_get_alpha
#' @return \code{sir_get_alpha_plusone(agent_states, model_config)} returns a vector, whose entries are alpha_i or alpha_r or 0 when previous states are 0, 1, 2 respectively.
#' @export
sir_get_alpha_plusone <- function(agent_states, model_config){
  alpha_plusone <- rep(0, length(agent_states));
  alpha_plusone <- (agent_states == 0) * (model_config$lambda * (agent_states == 1) %*% model_config$adjacency_matrix_b) + (agent_states == 1) * (model_config$gamma);
  # alpha_plusone <- (agent_states == 0) * pmax(1e-4, model_config$lambda * (agent_states == 1) %*% model_config$adjacency_matrix_b) + (agent_states == 1) * (model_config$gamma);
  return(as.vector(alpha_plusone))
}
