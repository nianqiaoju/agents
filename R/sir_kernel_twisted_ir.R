#' @title Twisted SIR kernels
#' @description  Function to sample from xt given it, rt and x(t-1), used in sir_csmc.R.
#' @param xxprev agent states at the previous time;
#' @param it infection counts at the current time;
#' @param rt recovery counts at the current time;
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @return a sample of xxprev, vector of length N
#' @export

sir_kernel_twisted_ir <- function(xxprev, it, rt, model_config){
  ## need to add some error messages to see if xxprev, it, rt are compatible
  alpha_matrix <- sir_get_alpha_sir(xxprev, model_config);
  xxnew <- rep(NA, N)
  recovery_prev <- which(xxprev == 2)
  xxnew[recovery_prev] <- 2 ## recovered stay there 
  rprev <- length(recovery_prev)
  ## only old infections can become new recovery
  previously_infected_agents <- which(xxprev == 1)
  iprev <- length(previously_infected_agents);
  s2i_count <- it - iprev + rt - rprev;
  ## check if conditions are met 
  if ( it + rt > model_config$N) warning('it + rt too big')
  if ( rt - rprev > iprev){
    warning('too many recovery')
  }
  if ( s2i_count > model_config$N - iprev - rprev) warning('too many new infections')
  if ( s2i_count < 0 ) warning('eta1 too small')
  newly_recovered_agents <- rcondbern(sum_x = rt - rprev, alpha = alpha_matrix[3,previously_infected_agents]); 
  ## newly_recovered is a binary vector of length iprev and 1 indicates transition from i to r 
  ## newly_recovered_agents goes from 1 to 2 and previous_infected_agents but not newly_recovered_agents stays at 1
  xxnew[previously_infected_agents] <- 1 + newly_recovered_agents
  ## only previously susceptible individuals can be infected
  previously_susceptible_agents <- which(xxprev == 0)
  infection_probability <- alpha_matrix[2, previously_susceptible_agents]; 
  newly_infected_agents <- rcondbern(sum_x = s2i_count, alpha = infection_probability);
  xxnew[previously_susceptible_agents] <- newly_infected_agents
  if (any(xxnew < xxprev)) print('twisted sir kernel error')
  if (sum(xxnew == 1) != it){
    warning('twisted sir kernel error 1')
    print(xxprev)
    print(it)
    print(rt)
    print(s2i_count)
    print(rt - rprev)
  }
  if (sum(xxnew == 2) != rt){
    warning('twisted sir kernel error 2')
    print(xprev)
    print(eta1)
    print(eta2)
    print(s2i_count)
    print(rt - rprev)
  } 
  return(xxnew)
}

