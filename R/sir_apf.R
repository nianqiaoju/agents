#' @title Auxilliary particle filter for the agent-based hidden SIR process.
#' @description unbiased approximation of the marginal likelihood. For this function, the observation density is \deqn{y_t \mid x_t \sim \textrm{Bin}(I(x_t),\rho).}
#' @param y a list of observations, length (T+1).
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @param particle_config a list specifying the particle filter, must pass the test of check_particle_config;
#' @return the particle filter outcome, must contain estimate of log marginal likelihood and effective sample size.
#' @export

sir_apf <- function(y, model_config, particle_config){
	num_observations <- length(y);
  loglikelihood <- rep(-Inf,num_observations);
  ess <- rep(1, num_observations);
  ancestors <- c(1:particle_config$num_particles);
  logweights <- logw <- rep(-Inf, particle_config$num_particles)
  logW <- rep(log(1 / particle_config$num_particles), particle_config$num_particles);
  ### logw: cumulative weights, used for resampling
  ### logW: log of normalized logw
  ### logweights: potential at the current step, use these to compute normalizing constant
  if(particle_config$save_particles) particles <- array(NA, dim = c(model_config$N, num_observations, particle_config$num_particles))
  if(particle_config$save_genealogy){
    eves <- c(1 : particle_config$num_particles);
  }
  if(particle_config$clock){
    runtimes <- rep(NA,num_observations);
    startstopwatch <- as.numeric(Sys.time());
  }
  # treat t = 0 differently 
  current_support_y <- y[0 + 1] : model_config$N;
  if (all(model_config$alpha0 == model_config$alpha0[1])){
    if(particle_config$verbose) cat("homogenous intial infection probability \n")
    logvt <- dbinom(x = current_support_y, size = model_config$N , prob = model_config$alpha0[1], log = TRUE) + 
      dbinom(x = y[0 + 1], size = current_support_y, prob = model_config$rho, log = T);
  }else{
    logvt <- logdpoisbinom_cpp(alpha = model_config$alpha0)[current_support_y + 1] + 
      dbinom(x = y[0 + 1], size = current_support_y, prob = model_config$rho, log = T);
  }
  loglikelihood[0 + 1] <- lw.logsum(logvt); ## p(y0)
  if(particle_config$clock) runtimes[1] <- Sys.time()[1];
  ## normalize to get conditional density of i0 | y0 
  vt <- lw.normalize(logvt);
  ## sample i0 | y0
  its <- sample(current_support_y, size = particle_config$num_particles, replace = T, prob = vt);
  ## sample x0 | i0
  if (all(model_config$alpha0 == model_config$alpha0[1])){
    if(particle_config$verbose) cat("homogenous intial infection probability \n");
    xts <- matrix(0, nrow = model_config$N, ncol = particle_config$num_particles); ## this line initializes the particles
    for (p in 1:particle_config$num_particles){
      xts[sample(x = model_config$N, size = its[p], replace = F) , p] <- 1;
    }
  }else{
    xts <- sapply(its, function(ii) rcondbern(sum_x = ii, alpha = model_config$alpha0));
  }
  if(particle_config$save_particles)  particles[,1,] <- xts;
  ## for t = 0, no need to resample because weights are the same
  for (t in 1:(num_observations - 1)){
    ## calculate p(x[t+1,k] = 1 | xt)
    alphats <- apply(xts, MARGIN = 2, FUN = function(xx) sir_get_alpha_i(xx,model_config))
    ## calculate the weights p(y[t+1] | xt) and conditional pmf p(i[t+1]| xt, y[t+1]) 
    current_support_y <- y[t + 1] : model_config$N
    logvt <- matrix(data = NA, nrow = length(current_support_y), ncol = particle_config$num_particles)
    logweights <- rep(NA, particle_config$num_particles)
    for (p in 1 : particle_config$num_particles){
      logvt[,p] <- logdpoisbinom_cpp(alpha = alphats[,p])[current_support_y + 1] + dbinom(x = y[t + 1], size = current_support_y, prob = model_config$rho, log = TRUE);
    }
    logweights <- lw_logsum_normalize_byCol(logvt);## logvt has been normalized in this step
    logweights[is.na(logweights)] <- - Inf
    loglikelihood[t + 1] <- lw.logsum(logweights) - log(particle_config$num_particles);
    if(particle_config$clock) runtimes[t+1] <- Sys.time()
    if(all(is.infinite(logweights))){## if particles degenerate
      break
    }
    ## resample: normalize the log weights
    weights <- lw.normalize(logweights)
    ess[t] <- 1 / sum(weights**2) / particle_config$num_particles;
    ## adaptive resampling
    if(ess[t] < particle_config$ess_threshold | any(weights == 0)){
      ancestors <- sample.int(n = particle_config$num_particles, prob = weights, replace = TRUE);
      logW <- rep(-log(particle_config$num_particles), particle_config$num_particles);
      logw <- rep(0, particle_config$num_particles);
      if(particle_config$verbose) cat("RESAMPLE!\ness at step t =", t, "is", ess[t],"\n");
    }else{
      ancestors <- 1 : particle_config$num_particles;
      logW <- log(weights);
      if(particle_config$verbose) cat("ess at step t =", t, "is", ess[t],"\n");
    }
    ## sample it given y and sample xt given it
    xnew <- matrix(NA, nrow = model_config$N, ncol = particle_config$num_particles)
    if (y[t + 1] == model_config$N){
      its <- rep(N, particle_config$num_particles)
      xnew <- matrix(1, nrow = model_config$N, ncol = particle_config$num_particles)
    }else{
      # sample its with ancestors
      its <- sapply(1:particle_config$num_particles, function(iparticle) sample(x = current_support_y, size = 1, prob = logvt[,ancestors[iparticle]]));
      xnew <- sapply(1:particle_config$num_particles, function(iparticle)
      	sir_kernel_twisted_i(xts[,ancestors[iparticle]], its[iparticle], model_config));
    }
    xts <- xnew
    if(particle_config$save_genealogy & particle_config$verbose) cat(ancestors, '\n')
    if(particle_config$save_genealogy) eves <- eves[ancestors]
    if(particle_config$save_particles){
      particles[,,] <- particles[, ,ancestors]
      particles[,t+1,] <- xts
    }
  }
  result <- list(log_incremental_likelihood = loglikelihood, 
                 log_marginal_likelihood = cumsum(loglikelihood), 
                 log_final_likelihood = sum(loglikelihood),
                 ess = ess);
  if(particle_config$clock){
    runtimes <- runtimes - startstopwatch
    result$runtimes = runtimes
    result$totaltime = runtimes[length(runtimes)]
  }
  if(particle_config$save_particles) result$particles = particles
  if(particle_config$save_genealogy) result$eves = eves
  return(result)
}
