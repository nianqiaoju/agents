#' @title SIR model likelihood
#' @description p(x,y) for the SIR model
#' @param xx agent states
#' @param y population observations 
#' @param model_config a list containinng parameters, features, and network structure
#' @export
#' @return log likelihood

sir_loglikelihood_complete <- function(xx, y, model_config){
	## first make sure xx and y are compatible
	if(dim(xx)[2] != length(y)) stop("agent states does not match observations");
	num_observations <- length(y);
	## logf0
	loglik <- sir_log_transition_density(xx[,1], model_config$alpha0);
	## log p(xt|x[t-1]):
	for(t in 1 : (num_observations -1)){
		loglik <- loglik + sir_log_transition_density(xx[,t + 1] - xx[,t], sir_get_alpha_plusone(xx[,t], model_config));
	}
	return(loglik);
}