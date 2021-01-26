#' @title SIR blocked gibbs 
#' @description update the agent states for the SIR model given observations such that the complete likelihood is not zero.
#' @param xx a matrix of size model_config$N by length(y), can be initialized from sis_init_hidden_popobservation.R 
#' @param block_members a list of agent ids, length should not exceed 5
#' @param time_step 0 <=t <= T
#' @param model_config a list that must contain
#' @param all_xx_block default is NULL. This matrix represents the set {1,2,3}^B.
#' \itemize{
#' \item N : population size
#' \item rho : reporting rate
#' \item alpha0 : vector of probabilities
#' \item state_to_alpha: function to compute alpha given agent states at a previous time 
#' }
#' @return xx

sir_xx_gibbs_blocked_kernel <- function(xx, block_members, y, model_config, all_xx_block = NULL){
	if(length(block_members)>5) stop("block size too large for exact updates");
	if(is.null(all_xx_block)) all_xx_block <- sir_get_state_space(length(block_members));
	## compute some quantities that are used repeatedly
	sum_xx_rest <- colSums(xx[-block_members,]==1);
	sum_all_xx_block <- rowSums(all_xx_block == 1);
	## forward algoritm
	m <- matrix(NA, nrow = dim(all_xx_block)[1], ncol = length(y));
	## each row is a state and each column is a time
	## m[,t+1] = p(xt, y(0:t)) for t = 0,...,T
	logf0 <- apply(all_xx_block, 1, function(xx) sir_log_transition_density(xx, model_config$alpha0));
	logg0 <- dbinom(x = y[1], size = sum_xx_rest[1] + sum_all_xx_block, prob = model_config$rho, log = TRUE);
	m[,1] <- logf0 + logg0;
	## need to store all the markov transition kernel densities
	logft <- array(NA, dim = c(dim(all_xx_block)[1], dim(all_xx_block)[1], length(y) - 1));
	## logft is the transition kernel from t - 1 to t
	for (t in 1 : (length(y) - 1)){
		xx_prev <- xx[, t - 1 + 1];
		xx_t <- xx[, t + 1];
		for (istate in 1 : dim(all_xx_block)[1]){
			xx_prev[block_members] <- all_xx_block[istate,];
			alpha_prev2t <- sir_get_alpha_plusone(xx_prev, model_config);
			for(jstate in 1 : dim(all_xx_block)[1]){
				xx_t[block_members] <- all_xx_block[jstate,];
				logft[istate, jstate, t] <- sir_log_transition_density(xx_t - xx_prev, alpha_prev2t);
			}
		}
		m[, t + 1] <- dbinom(y[t+1],sum_xx_rest[t+1] + sum_all_xx_block, model_config$rho, TRUE) + apply(logft[,,t] , 2, function(log_transition) lw.logsum(log_transition + m[,t]));
	}
	filters <- apply(m, 2, lw.normalize); 
	## each row is a states
	## filters[,t - 1] = log(xt | y_{0:t}) for t = 1,2,...,T
	## backward sampling starts from terminal time
	## sample the index in all_xx_block first
	bs_indices <- rep(NA, length(y))
	bs_indices[length(y)] <- sample.int(dim(all_xx_block)[1], 1, prob = filters[,length(y)]);
	for (t in (length(y) - 1) : 1) {
		logprob <- logft[,bs_indices[t + 1],t] + log(filters[, t]);
		bs_indices[t] <- sample.int(dim(all_xx_block)[1], 1, prob = lw.normalize(logprob));
	}
	## convert to the states
	for(t in 1 : length(y)){
		xx[block_members, t] <- all_xx_block[bs_indices[t],];
	}
	return(list(xx = xx, m = m)) ## output m for debug, can remove it in production code
}