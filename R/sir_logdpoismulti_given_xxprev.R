#' @title SIR model
#' @description the log density of \eqn{p(st,it | x[t-1])}. 
#' NOTE: if network is full, then we can use sir_csmc_update_f_matrix instead.
#' @param xxprev \eqn{x[t-1]}
#' @param model_config a list containing parameters and network structure 
#' @export

sir_logdpoismulti_given_xxprev <- function(xxprev, model_config){
  ## compute the transition probability from xx[t-1] to the aggregated state (st, it)
  qq <- matrix(-Inf, nrow = model_config$N + 1, ncol = model_config$N + 1);
  prev_infected_agents <- which(xxprev == 1);
  prev_susceptible_agents <- which(xxprev == 0);
  icount_prev <- length(prev_infected_agents);
  scount_prev <-length(prev_susceptible_agents);
  alphai <- sir_get_alpha_i(xxprev, model_config);
  di2i <- logdpoisbinom_cpp(alphai[prev_infected_agents]); ## this vector starts at 0
  ds2i <- logdpoisbinom_cpp(alphai[prev_susceptible_agents]);
  for (s2i in  0 : icount_prev){
    for (i2r in 0 : scount_prev){
      it <- s2i + i2i;
      st <- scount_prev - s2i;
      qq[st + 1, it + 1] <- di2r[i2i + 1] + ds2i[s2i + 1];
    }
  }
  return(qq);
}
