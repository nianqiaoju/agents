#' @title Gibbs sampler for SIR model
#' @description initialize agent states for SIR model given y such that the complete likelihood is not zero.
#' @param y population observations 
#' @param model_config a list containinng parameters, features, and network structure
#' @export
#' @return agent_states

sir_xx_initialize <- function(y,model_config){
  num_observations <- length(y);
  xx  <- matrix(0, nrow = model_config$N, ncol = num_observations);
  ## sample it - yt for every time
  it <- rpois(n = num_observations, lambda =  0.5 * model_config$N  * (1 - model_config$rho)) + y;
  ## make sure 1 <= it <= yt 
  it <- pmin(it, N); 
  it <- pmax(it, 1);
  ## make sure there are it infections at every t
  for (t in 1:num_observations){
    xx[sample.int(n = N, size = it[t]) , t] <- 1;
  }    
  ## fill the agent states such that it looks like 0000111112222
  for (n in 1 : model_config$N){
  	infected_times <- which(xx[n,] == 1);
  	if (length(infected_times)){
  		start <- min(infected_times);
  		end <- max(infected_times);
  		xx[n,c(start:end)]<- 1 ;
  		if (end < num_observations)  xx[n,(end + 1):(num_observations)]<- 2;
  	}
  }
  return(xx)
}