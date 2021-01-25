#' @title SIS blocked gibbs 
#' @description update the agent states for the SIS model given observations such that the complete likelihood is not zero.
#' @param xx a matrix of size model_config$N by length(y), can be initialized from sis_init_hidden_popobservation.R 
#' @param block_members a list of agent ids, length should not exceed 5
#' @param model_config a list that must contain
#' @param y observations
#' @param all_xx_block default is NULL, this is a matrix containing the state space {0,1}^B.
#' @return xx


sis_xx_gibbs_block_by_agents_kernel <- function(xx, block_members, y, model_config, all_xx_block = NULL){
	if(length(block_members) > 5) stop("block size too large for exact updates");
	if(is.null(all_xx_block)) all_xx_block <- sis_get_state_space(length(block_members));
	## compute some quantities that are used repeatedly
	sum_xx_rest <- colSums(xx[-block_members, ]);
	sum_all_xx_block <- rowSums(all_xx_block);
	## forward algorithm
	m <- matrix(NA, nrow = dim(all_xx_block)[1], ncol = length(y)); 
	## each row is a state
	## each column is a time
	## m[ , t + 1] = p(xt, y(0:t)) for t = 0, ..., T
	logf0 <- apply(all_xx_block, 1, function(block) logdbern_sum_cpp(model_config$alpha0[block_members], block));
	logg0 <- dbinom(x = y[1], size = sum_xx_rest[1] + sum_all_xx_block, prob = model_config$rho, log = TRUE);
	m[,1] <- logf0 + logg0
	## need to store all the markov transition kernel densities
	logft <- array(NA, dim = c(dim(all_xx_block)[1], dim(all_xx_block)[1], length(y) - 1))
	## logft[,,t] is the transition kernel from t - 1 to t
	for (t in 1 : (length(y) - 1)){
		xx_prev <- xx[,t - 1 + 1];
		xx_t <- xx[, t + 1];
		for (istate in 1 : dim(all_xx_block)[1]){
			xx_prev[block_members] <- all_xx_block[istate,];
			alpha_prev2t <- sis_get_alpha(xx_prev, model_config);
			for (jstate in 1 : dim(all_xx_block)[1]){
				xx_t[block_members] <- all_xx_block[jstate, ];
				logft[istate,jstate,t] <- logdbern_sum_cpp(alpha_prev2t, xx_t);
			}
		}
		m[,t + 1] <- dbinom(y[t + 1], sum_xx_rest[t + 1] + sum_all_xx_block, model_config$rho, TRUE) + apply(logft[,,t], 2, function(log_transition) lw.logsum(log_transition + m[, t]));
	}
	filters <- apply(m, 2, lw.normalize);
	## each row is a state and filters[, t - 1]  = log p (xt | y_{0:t})
	## backward sampling 
	## sample the indices first 
	bs_indices <- rep(NA, length(y))
	bs_indices[length(y)] <- sample.int(dim(all_xx_block)[1], 1, prob = filters[,length(y)]);
	for (t in (length(y) - 1) : 1) {
		logprob <- logft[,bs_indices[t + 1],t] + log(filters[, t]);
		bs_indices[t] <- sample.int(n = dim(all_xx_block)[1], size = 1, replace = FALSE, prob = lw.normalize(logprob));
	}
	## convert incides to samples
	for(t in 1 : length(y)){
		xx[block_members, t] <- all_xx_block[bs_indices[t],];
	}
	return(list(xx = xx, m = m))
}
