#' @title SIR model transition kernels
#' @rdname sir_get_alpha
#' @description The alpha function for SIR process gives the transition probabilities from \eqn{x[t,n]} to \eqn{x[t+1,n]} given \eqn{x[t]}.
#' @param agent_states a vector containing agent states at time \eqn{t}.
#' @param model_config a list containing model parameters, which must pass the function \code{check_model_config()}.
#' @return \code{sir_get_alpha_sir(agent_states,model_config)} returns a matrix with three rows : alpha_s, alpha_i and alpha_r.
#' @export
#' 
sir_get_alpha_sir <- function(agent_states, model_config){
	alpha_matrix <- matrix(NA, nrow = 3, ncol = model_config$N);
	## probability of recovery
	alpha_matrix[3,] <- (agent_states == 2) + (agent_states == 1) * model_config$gamma
	alpha_matrix[2,] <- sir_get_alpha_i(agent_states, model_config);
	alpha_matrix[1,] <- 1 - alpha_matrix[3,] - alpha_matrix[2,];
	return(alpha_matrix)
}

#' @rdname sir_get_alpha
#' @return \code{sir_get_alpha_i(agent_states, model_config)} returns alpha_i as a vector.
#' @export
sir_get_alpha_i <- function(agent_states, model_config){
  if(model_config$network_type == "full"){
    it <- sum(agent_states == 1);
    alpha <- (agent_states == 0) *  (model_config$lambda * it / N) + (agent_states == 1) * (1 - model_config$gamma);
  }else{
    alpha <- (agent_state == 1) * ( 1 - model_config$gamma) + (agent_state == 0) * (model_config$lambda * (agent_state %*% model_config$adjacency_matrix_b));
  }
	return(as.vector(alpha));
}


#' @rdname  sir_get_alpha
#' @return \code{sir_get_alpha_plusone(agent_states, model_config)} returns a vector, whose entries are the probabilities that p(x[t+1,k] = x[t,k] + 1);
#' @export
sir_get_alpha_plusone <- function(agent_states, model_config){
  if(model_config$network_type == "full"){
    it <- sum(agent_states == 1);
    alpha_plusone <- (agent_states == 0) *  (model_config$lambda * it / N) + (agent_states == 1) * (model_config$gamma);
  }else{
    alpha_plusone <- (agent_state == 1) * (model_config$gamma) + (agent_state == 0) * (model_config$lambda * (agent_state %*% model_config$adjacency_matrix_b));
  }
  return(as.vector(alpha_plusone));
}
