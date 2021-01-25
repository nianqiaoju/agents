#' @title SIS swap kernel
#' @description update the agent states for the SIS model given observations such that the complete likelihood is not zero.
#' @param xx a matrix of size model_config$N by length(y), can be initialized from sis_init_hidden_popobservation.R 
#' @param time_step 0 <=t <= T
#' @param model_config a list that must contain
#' \itemize{
#' \item N : population size
#' \item rho : reporting rate
#' \item alpha0 : vector of probabilities
#' \item state_to_alpha: function to compute alpha given agent states at a previous time 
#' }
#' @return xx


sis_xx_gibbs_swap_at_time_kernel <- function(xx, time_step, model_config){
  if((time_step + 1) > dim(xx)[2]) stop("time exceeds observation length")
	xx_current <- xx[,time_step + 1]
	if (all(xx_current == 0) | all(xx_current == 1)) return(xx); ## no swap necessary if all 0s or 1s
	## choose indices for swap
	## NOTE: this implementation is O(N)
	## NOTE: swap moves do not change sum(xt) so the density does not depend on y[t]
	i0 <- sample(which(xx_current == 0), size = 1);
	i1 <- sample(which(xx_current == 0), size = 1);
	## previous step to now 
	if (time_step == 0){
		alpha_prev2now <- model_config$alpha0;
	}else{
		alpha_prev2now <- sis_get_alpha(xx[,time_step - 1 + 1], model_config);
	}
	log_accept_prob <- log(alpha_prev2now[i0]) + log(1 - alpha_prev2now[i1]) - log(alpha_prev2now[i1]) - log(1 - alpha_prev2now[i0]);
	## now to next step
	if ((time_step + 1) != dim(xx)[2]){
		xx_swap <- xx_current 
		xx_swap[i1] <- 0 
		xx_swap[i0] <- 1
		alpha_now2next_noswap <- sis_get_alpha(xx_current, model_config);
		alpha_now2next_swap <- sis_get_alpha(xx_swap, model_config);
		log_accept_prob <- log_accept_prob + logdbern_sum_cpp(alpha_now2next_swap, xx[, time_step + 1]) - logdbern_sum_cpp(alpha_now2next_noswap, xx[, time_step + 1]);
	}
	## accept step
	if (log(runif(1)) < log_accept_prob){
		## swap 
		xx[i0, time_step + 1] <- 1
		xx[i1, time_step + 1] <- 0 
	}
	return(xx)
}