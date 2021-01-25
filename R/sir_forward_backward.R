#' Forward backward algorithm for the SIR model
#' @param y population observations
#' @param model_config a list containing
#' \itemize{
#' \item 'N' : population size
#' \item 'alpha0'
#' \item 'lambda'
#' \item 'gamma'
#' \item 'rho'
#' \item 'adjacency_matrix_b': network structure
#' }
#' @param num_samples number of backward samples needed, set to 1 if not specifed. If num.samples = 0, then no backward sample is computed.
#' @param all_xx a 2**N by N matrix, representing the sample space, pre-calculated to save time 
#' @return a list containing 
#' \itemize{
#' \item 'log marginal likelihood and one backward sample' 
#' \item 'bs' : backward samples in a [N * num_observations * num_samples] array
#' }
#' @export

sir_forward_backward <- function(y, model_config, num_samples = 1, all_xx = NULL){
	if(model_config$N > 5) stop('population too large for forward backward');
	if(is.null(all_xx)) all_xx <- sir_get_state_space(model_config$N);
	num_observations <- length(y);
	## f0 density
	logf0 <- rep(NA, dim(all_xx)[1])
	logf0 <- apply(all_xx, 1, function(xx) sir_log_transition_density(xx, model_config$alpha0));
	## log f(x(t) | x(t-1))
	all_alpha_plusone <- t(apply(all_xx, 1, function(xx) sir_get_alpha_plusone(xx, model_config)))
	## log markov transition kernel between states
	logdtransition <- matrix(NA, nrow = dim(all_xx)[1], ncol = dim(all_xx)[1])
	for(istate in 1 : dim(all_xx)[1]){
		logdtransition[istate,] <- apply(all_xx, 1, function(xx) sir_log_transition_density(xx - all_xx[istate,], all_alpha_plusone[istate,]));
	}
	# apply(logdtransition, 1, lw.logsum); ## DEBUG CODE
	## I(xx) for each state
	all_sum_x <- rowSums(all_xx == 1)
	# each row correspond to one individual
	# each column correspond to one day 
	# M[,1] = log p(x0,y0) 
	# M[,t + 1] = log p(x_t,y_{0:t}) for all other columns // these are .R indices
	M <- sis_forward_algorithm_cpp(logf0, logdtransition, y, all_sum_x, model_config$rho);
	log_marginal_likelihood <- apply(M, 2, lw.logsum); ## log(p(y_{0:t}))
	log_incremental_likelihood <- rep(NA, num_observations);
	log_incremental_likelihood[1] <- log_marginal_likelihood[1];
	log_incremental_likelihood[-1] <- diff(log_marginal_likelihood, lag = 1);
	result <- list(M = M, log_marginal_likelihood = log_marginal_likelihood, log_incremental_likelihood = log_incremental_likelihood)
	## backward sampling steps::
	if (num_samples >= 1){
		## sample the index of states in the all_xx matrix
		## each row represents one posterior sample
		bs_indices <- matrix(NA, nrow = num_samples, ncol = num_observations); 
		## start at the terminal time
		bs_indices[, num_observations] <- sample.int(n = dim(all_xx)[1], size = num_samples, replace = TRUE, prob = lw.normalize(M[, num_observations]));
		## start the backward recursion
		for( t in (num_observations - 2) : 0){
			bs_indices[, t + 1] <- sapply(bs_indices[, t + 1 + 1], 
				function(istate) sample.int(n = dim(all_xx)[1], size = 1, prob = lw.normalize(logdtransition[, istate] + M[, t + 1 + 1])));
		}
		## convert the indices to binary vectors
		bs <- array(dim = c(num_samples, dim(all_xx)[2], num_observations));
		for(isample in 1 : num_samples){
			for(t in 1 : num_observations){
				bs[isample, , t] <- all_xx[bs_indices[isample, t], ]
			}
		}
		result$bs = bs
	}
	return(result)
}