#' @title Markov transition kernel for SIR model
#' @description computes the log transition density f(x[t+1] | x[t])
#' @param xx_diff binary vector, equals x[t+1] - x[t]
#' @param alpha_plueone probability of transitioning to the next state from x[t]
#' @export
sir_log_transition_density <- function(xx_diff, alpha_plusone){
	if(any(xx_diff != 1 & xx_diff != 0) ) return(-Inf);
	return(logdbern_sum_cpp(alpha_plusone, xx_diff));
}