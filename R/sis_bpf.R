#' @title Bootstrap Particle Filter for SIS model with Population Observations
#' @description approximates marginal likelihood of the SIS model using bootstrap particle filter
#' NOTE: this cannot do predictions yet.
#' @param y a vector of length (T+1)
#' @param model_config a list containing:
#' \itemize{
#' \item 'N': size of the population
#' \item 'alpha0' : initial infection probability
#' \item 'lambda': infection rate 
#' \item 'gamma': recovery rate
#' \item 'adjacency_matrix_b': network structure
#' \item 'rho': reporting rate
#' \item }
#' @param particle_config a list containing:
#' \itemize{
#' \item 'num_particles' : number of particles for the bootstrap particle filter
#' \item 'ess_threshold' : if effective sample size drops below the threshold, then perform a resample step. ess_threshold = 1 means resampling at every step.
#' \item 'save_particles': binary
#' \item 'clock' : binary, default to FALSE. If clock = TRUE, then we will use a stopwatch to document its Sys.time()
#' \item 'save_genealogy': binary
#' \item 'verbose': if TRUE print messages for debug
#' }
#' @return A list containing 
#' \itemize{
#' \item 'log_incremental_likelihood' :  
#' \item 'log_marginal_likelihood' :
#' \item 'log_final_likelihood' :
#' \item 'particles': [N, (T+1), num_particles] array storing the particles: 
#' \item 'eves' : the unique ancesters of each particle in the genealogy structure, if save_genealogy = TRUE;
#' \item 'ess' : effective sample size at each step
#' \item 'runtime' : elapse time at each step, if clock = TRUE;
#' \item 'totaltime' : elapse time of running the whole particle filter, if clock = TRUE;
#' }
#' @export


sis_bpf <- function(y,model_config,particle_config){
  num_observations <- length(y)
  loglikelihood <- rep(-Inf,num_observations)
  if(particle_config$save_particles)  particles <- array(NA, dim = c(model_config$N, num_observations, particle_config$num_particles));
  if(particle_config$save_genealogy){
    eves <- c(1 : particle_config$num_particles)
  } 
  if(is.null(particle_config$ess_threshold)) stop("Please specify ess_threshold for resampling.");
  if(particle_config$verbose) cat("[ess threhold is", particle_config$ess_threshold,"]\n");
  logweights <- logw <- rep(0, particle_config$num_particles)
  logW <- rep(log(1 / particle_config$num_particles), particle_config$num_particles);
  ### logw: cumulative weights, used for resampling
  ### logW: log of normalized logw
  ### logweights: potential at the current step, use these to compute normalizing constant
  ess <- rep(0, num_observations)
  ancestors <- c(1:particle_config$num_particles)
  if(particle_config$clock){
    runtimes <- rep(NA,num_observations)
    startstopwatch <- as.numeric(Sys.time())
  }
  ## start from t = 0
  alphats <- matrix(model_config$alpha0, nrow = model_config$N, ncol = particle_config$num_particles);
  for (t in 1 : num_observations){
    ## sample according transition from the model for each particle
    xts <- (runif(particle_config$num_particles * model_config$N) < alphats);
    ## compute weights by p(yt | xt)
    its <- colSums(xts);
    logweights <- dbinom(x = as.numeric(y[t]), size = its, prob = model_config$rho, log = TRUE);
    maxlogweights <- max(logweights);
    ## if the particle filter collapses we should exit early and print an error message (?)
    if(is.infinite(maxlogweights)){
      warning("BPF for SIS collapses");
      break
    }
    logweights <- logweights - maxlogweights;
    loglikelihood[t] <- log(mean(exp(logweights))) + maxlogweights;
    if(particle_config$clock) runtimes[t] <- as.numeric(Sys.time());
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    ess[t] <- 1 / sum(weights ** 2) / particle_config$num_particles;
    ## resample steps 
    if(ess[t] < particle_config$ess_threshold){
      ancestors <- sample.int(n = particle_config$num_particles, prob = weights, replace = TRUE);
      logW <- rep(-log(particle_config$num_particles), particle_config$num_particles);
      logw <- rep(0, particle_config$num_particles);
      if(particle_config$verbose) cat("RESAMPLE!\ness at step t =", t, "is", ess[t],"\n");
    }else{
      ancestors <- 1 : particle_config$num_particles;
      logW <- log(weights);
      if(particle_config$verbose) cat("ess at step t =", t, "is", ess[t],"\n");
    }
    xts[,] <- xts[,ancestors];
    ## store genealogy
    if(particle_config$save_genealogy) eves <- eves[ancestors];
    ## store particles
    if(particle_config$save_particles){
      particles[,,] <- particles[, ,ancestors]
      particles[,t+1,] <- xts
    }
    ## update alpha for next step
    alphats <- apply(xts, MARGIN = 2, FUN = function(xx) sis_get_alpha(agent_state = xx, model_config = model_config)); ## p(X(t+1,k)=1|xt) 
  }
  result <- list(log_incremental_likelihood = loglikelihood, 
                 log_marginal_likelihood = cumsum(loglikelihood), 
                 log_final_likelihood = sum(loglikelihood),
                 ess = ess)
  if(particle_config$clock){
    runtimes <- as.numeric(runtimes - startstopwatch, units = "seconds");
    result$runtimes = runtimes
    result$totaltime = runtimes[length(runtimes)]
  }
  if(particle_config$save_particles) result$particles = particles
  if(particle_config$save_genealogy) result$eves = eves
  return(result)
}
