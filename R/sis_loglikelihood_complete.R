#' @title Complete Likelihood for SIS model
#' @description computes \eqn{log p(x, y)} when agent states and counts are observed
#' @param y a vector of observations
#' @param agent_state binary matrix of size N by length(y)
#' @param model_config a list containing parametesr and network structure, passing the check of \code{check_model_config()}.
#' @export
#' 
sis_loglikelihood_complete <- function(y, agent_states, model_config){
  num_observations <- length(y);
	## check if observations and hidden states are compatible
	if (num_observations != dim(agent_states)[2]) warning('incorrect length of observations');
	if (N != dim(agent_states)[1]) warning('incorrect number of agents');
	## observation densities
	llik <- sum(dbinom(x = y, size = colSums(agent_states), prob = model_config$rho, log = TRUE ));
	## transition densities
	llik <- llik + logdbern_sum_cpp(model_config$alpha0, agent_states[,1]); ## t = 0
	for (t in 1 : (num_observations - 1)){## t = 1,...,T
		a <- sis_get_alpha(agent_states[, t - 1 + 1], model_config);
		llik <- llik + logdbern_sum_cpp(a, agent_states[,t + 1]);
	}
	return(llik);
}