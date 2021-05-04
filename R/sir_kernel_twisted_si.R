#' @title Twisted SIR kernels 
#' @description  Function to sample from xt given it, st and x(t-1), used in sir_csmc.R.
#' @param xxprev agent states at the previous time;
#' @param  st infection counts at the current time;
#' @param it infection counts at the current time;
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @return a sample of xxprev, vector of length N
#' @export
#' 
#' 
sir_kernel_twisted_si <- function(xxprev, st, it, model_config){
  alpha_i <- sir_get_alpha_i(xxprev, model_config); # individual transition probabilities
  xxnew <- integer(model_config$N); # stores the next states
  recovery_prev <- which(xxprev == 2);
  rcount_prev <- length(recovery_prev);
  xxnew[recovery_prev] <- 2; ## recovered stay there 
  ## only old infections can become new recovery
  previously_infected_agents <- which(xxprev == 1);
  icount_prev <- length(previously_infected_agents);
  scount_prev <- model_config$N - icount_prev - rcount_prev;
  s2i_count <- scount_prev - st;
  ## check if conditions are met 
  if (it + st > model_config$N) warning('it + st too big');
  if (scount_prev + icount_prev - st - it > icount_prev) warning('too many recovery');
  if (s2i_count > scount_prev - st) warning('too many new infections');
  if (s2i_count < 0) warning('not enough new infections');
  ## newly_recovered is a binary vector of length iprev and 1 indicates transition from i to r 
  ## newly_recovered_agents goes from 1 to 2 and previous_infected_agents but not newly_recovered_agents stays at 1
  newly_recovered_agents <- rcondbern(sum_x = scount_prev + icount_prev - st - it, alpha = 1 - alpha_i[previously_infected_agents]);
  xxnew[previously_infected_agents] <- 1 + newly_recovered_agents;
  ## only previously susceptible individuals can be infected
  previously_susceptible_agents <- which(xxprev == 0);
  infection_probability <- 1 - alpha_i[previously_susceptible_agents]; 
  newly_infected_agents <- rcondbern(sum_x = s2i_count, alpha = infection_probability);
  xxnew[previously_susceptible_agents] <- newly_infected_agents
  if (any(xxnew < xxprev)) print('twisted sir kernel error')
  if (sum(xxnew == 1) != it){
    warning('twisted sir kernel error 1')
    print(xxprev)
    print(it)
    print(rt)
    print(s2i_count)
  }
  if (sum(xxnew == 0) != st){
    warning('twisted sir kernel error 2')
    print(xxprev)
    print(eta1)
    print(eta2)
    print(s2i_count)
  } 
  return(xxnew)
}
