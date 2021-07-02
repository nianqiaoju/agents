sis_loglikelihood_marginal <- function(y, model_config, particle_config, method){
	method <- match.arg(method, choices = c('apf','exact','csmc'))
	## update model_config according to parameters
	if(method == 'mc' & is.null(particle_config$num_particles)){
    stop("please specify the number of particles.")
  }
  if(method == 'exact' & model_config$N >= 5) stop("population too big for exact computations");
  if (method == 'exact'){
  	loglik <- sis_forward_backward(y, model_config, 0, state_space)$log_marginal_likelihood[length(y)]
  	return(loglik)
  }
  if (method == 'csmc'){
    loglik <- sis_csmc(y, model_config, particle_config)$log_final_likelihood
    return(loglik)
  }
  if(method == 'apf'){
  	loglik <- sis_apf(y, model_config, particle_config)$log_final_likelihood
  	return(loglik)
  }
	return(-Inf)
}