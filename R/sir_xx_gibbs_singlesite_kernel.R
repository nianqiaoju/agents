#' @title Gibbs sampler for SIR model
#' @param xx current agent states
#' @param y observations
#' @param model_config a list containing paraemters, features and network structure
#' @return xx
#' @export 

sir_xx_gibbs_singlesite_kernel <- function(xx, y, model_config){
	num_observations <- length(y);
	for(n in 1 : model_config$N){
		infected_times <- which(xx[n,] == 1);
		if (length(infected_times) == 0){
			## if this agent is never infected, can only be infected at terminal time
			xx0 <- xx1 <- xx[,num_observations];
			xx1[n] <- 1;
			a <- sir_get_alpha_i(xx[,num_observations - 1], model_config)[n];
			sum_xx0 <- sum(xx0 == 1);
			lp0 <- log(1 - a) + dbinom(y[num_observations], sum_xx0, model_config$rho, TRUE);
			lp1 <- log(a) + dbinom(y[num_observations], sum_xx0 + 1, model_config$rho, TRUE);
			maxlp <- max(lp0, lp1);
			lp0 <- lp0 - maxlp;
			lp1 <- lp1 - maxlp;
			if (runif(1) < exp(lp1) / (exp(lp1) + exp(lp0))){
				xx[n, num_observations] <- 1;
				start <- num_observations;
			}
		}else{
			start <- infected_times[1];
			end <- mylast(infected_times);
			## first consider the s2i transition
			## xx[n, start - 1] = 0 or 1
			if(start >= 2){
				if(start == 2){ 
					a <- model_config$alpha0[n];
				}
				if(start > 2){
					a <- sir_get_alpha_i(xx[,start - 2], model_config)[n];
				}
				xx0 <- xx[,start - 1]; sum_xx0 <- sum(xx0 == 1);
				xx1 <- xx[,start - 1]; xx1[n] <- 1;
				lp0 <- log(1 - a) +  dbinom(y[start - 1],sum_xx0, model_config$rho, TRUE) + sir_log_transition_density(xx[,start] - xx0, sir_get_alpha_plusone(xx0,  model_config));
				lp1 <- log(a) + dbinom(y[start - 1], sum_xx0 + 1, model_config$rho, TRUE) + sir_log_transition_density(xx[,start] - xx1, sir_get_alpha_plusone(xx1, model_config)) ;
				maxlp <- max(lp0, lp1);
				lp0 <- lp0 - maxlp;
				lp1 <- lp1 - maxlp;
				if (runif(1) < exp(lp1) / (exp(lp1) + exp(lp0))){
					xx[n, start - 1] <- 1
					start <- start - 1
				}
			}
			## xx[n,start] = 0 or 1 and if xx[n,start] <-0, start++, need to consider the next one
			finish_start_scan <- FALSE;
			while(!finish_start_scan){
				if(start != end){
					xx0 <- xx1 <- xx[,start];
					xx0[n] <- 0;
					sum_xx1 <- sum(xx1 == 1);
					if(start == 1){
						a <- model_config$alpha0[n];
					}else{
						a <- sir_get_alpha_i(xx[,start - 1], model_config)[n];
					}
					lp0 <- log(1 - a) + dbinom(y[start],sum_xx1 - 1,model_config$rho, TRUE);
					lp1 <- log(a) +  dbinom(y[start],sum_xx1,model_config$rho, TRUE);
					if (start < num_observations){
						lp0 <- lp0 + sir_log_transition_density(xx[,start + 1] - xx0, sir_get_alpha_plusone(xx0, model_config));
						lp1 <- lp1 + sir_log_transition_density(xx[,start + 1] - xx1, sir_get_alpha_plusone(xx1, model_config));
					}
					maxlp <- max(lp0, lp1);
					lp0 <- lp0 - maxlp;
					lp1 <- lp1 - maxlp;
					if (runif(1) < exp(lp1) / (exp(lp0) + exp(lp1))){
						finish_start_scan <- TRUE
					}else{
						xx[n, start] <- 0
						start <- start + 1;
					}
				}else{
					finish_start_scan <- TRUE;
				}
			}
			## consider the i2r transition
			finish_end_scan <- FALSE;
			if (start != end){
				## xx[n, end] = 1 or 2
				xx1 <- xx2<- xx[,end];
				xx2[n] <- 2;
				sum_xx1 <- sum(xx1 == 1);
				lp1 <- log(1 - model_config$gamma[n]) + dbinom(y[end],sum_xx1,model_config$rho, TRUE);
				lp2 <- log(model_config$gamma[n]) + dbinom(y[end],sum_xx1 - 1,model_config$rho, TRUE);
				if (end < num_observations){
					lp1 <- lp1 + sir_log_transition_density(xx[, end + 1] - xx1, sir_get_alpha_plusone(xx1, model_config));
					lp2 <- lp2 + sir_log_transition_density(xx[, end + 1] - xx2, sir_get_alpha_plusone(xx2, model_config));
				}
				maxlp <- max(lp1, lp2);
				lp1 <- lp1 - maxlp;
				lp2 <- lp2 - maxlp;
				if (runif(1) < exp(lp2) / (exp(lp1) + exp(lp2))){
					xx[n, end] <- 2;
					end <- end - 1;
					finish_end_scan <- TRUE;
				}
			}
			while(!finish_end_scan){
				## xx[n,end + 1] = 1 or 2;
				if (end + 1 <= num_observations){
					xx1 <- xx2 <- xx[,end + 1];
					xx1[n] <- 1;
					sum_xx1 <- sum(xx1 == 1);
					lp1 <- log(1 - model_config$gamma[n]) + dbinom(y[end + 1],sum_xx1,model_config$rho, TRUE);
					lp2 <- log(model_config$gamma[n]) + dbinom(y[end + 1],sum_xx1 - 1,model_config$rho, TRUE);
						if (end + 2 <= num_observations){
							lp1 <- lp1 + sir_log_transition_density(xx[, end + 2] - xx1, sir_get_alpha_plusone(xx1, model_config));
							lp2 <- lp2 + sir_log_transition_density(xx[, end + 2] - xx2, sir_get_alpha_plusone(xx2, model_config));
						}
					maxlp <- max(lp1, lp2);
					lp1 <- lp1 - maxlp;
					lp2 <- lp2 - maxlp;
					if (runif(1) < exp(lp2) / (exp(lp1) + exp(lp2))){
						xx[n, end + 1] <- 2;
						finish_end_scan <- TRUE;
					}
				}else{
					finish_end_scan <- TRUE;
				}
			}
		}
	}
	return(xx)
}