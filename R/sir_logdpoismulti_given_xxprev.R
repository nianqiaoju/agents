#' @title SIR model
#' @description the density of p(it, rt | x[t-1]) 
#' @param xxprev x[t-1] 
#' @param model_config a list containing parameters and network structure 
#' @export
#' 
sir_logdpoismulti_given_xxprev <- function(xxprev, model_config){
  ## compute the transition probability from xx[t-1] to (it,rt)
  qq <- matrix(NA, nrow = model_config$N + 1, ncol = model_config$N + 1);
  ## compute alpha_s and alpha_r and alpha_i for each agent
  alpha_matrix <- sir_get_alpha_sir(xxprev, model_config);
  ## p(rt | xx[t-1])
  prev_infected_agents <- which(xxprev == 1);
  prev_susceptible_agents <- which(xxprev == 0);
  prev_infection_counts <- length(prev_infected_agents);
  prev_recovery_counts <- model_config$N - length(prev_susceptible_agents) - prev_infection_counts;
  support_i2r_counts <- 0 : length(prev_infected_agents);
  support_s2i_counts <- 0 : length(prev_susceptible_agents);
  di2r <- logdpoisbinom_cpp(alpha_matrix[3,prev_infected_agents]);
  ds2i <- logdpoisbinom_cpp(alpha_matrix[2,prev_susceptible_agents]);
  for (s2i in support_s2i_counts){
    for (i2r in support_i2r_counts){
        it <- prev_infection_counts + s2i - i2r ;
        rt <- prev_recovery_counts + i2r;
        qq[it + 1, rt +1] <- di2r[i2r + 1] + ds2i[s2i + 1];
    }
  }
  qq[which(is.na(qq))] <- -Inf
  return(qq)
}


sir_logdpoismulti_given_xxprev_si <- function(xxprev, model_config){
  ## compute the transition probability from xx[t-1] to the aggregated state (st, it)
  qq <- matrix(-Inf, nrow = model_config$N + 1, ncol = model_config$N + 1);
  prev_infected_agents <- which(xxprev == 1);
  prev_susceptible_agents <- which(xxprev == 0);
  icount_prev <- length(prev_infected_agents);
  scount_prev <-length(prev_susceptible_agents);
  support_i2r_counts <- 0 : icount_prev;
  support_s2i_counts <- 0 : scount_prev;
  alphai <- sir_get_alpha_i(xxprev, model_config);
  di2r <- logdpoisbinom_cpp(1 - alphai[prev_infected_agents]); ## this vector starts at 0
  ds2i <- logdpoisbinom_cpp(alphai[prev_susceptible_agents]);
  for (s2i in support_s2i_counts){
    for (i2r in support_i2r_counts){
      it <- icount_prev + s2i - i2r ;
      st <- scount_prev - s2i;
      qq[st + 1, it + 1] <- di2r[i2r + 1] + ds2i[s2i + 1];
    }
  }
  return(qq);
}
