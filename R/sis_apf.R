#'@title Fully adapted auxiliary Particle Filter for SIS model with Population Observation
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
#' \item 'verbose':
#' \item 'exact': binary, if TRUE, then uses exact PoisBinom and exact CondBern, if FALSE then use translated Poisson and MCMC for CondBern. MUST HAVE!
#' \item 'num_mcmc": required if exact = FALSE;
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

sis_apf <- function(y, model_config, particle_config){
  num_observations <- length(y)
  loglikelihood <- rep(-Inf,num_observations)
  if(particle_config$save_particles)  particles <- array(NA, dim = c(model_config$N, num_observations, particle_config$num_particles))
  if(particle_config$save_genealogy){
    eves <- c(1 : particle_config$num_particles)
  } 
  if(is.null(particle_config$ess_threshold)) stop("Please specify ess_threshold for resampling.\n");
  if(is.null(particle_config$exact)) stop("Please specify method to for density evaluations and for sampling.\n")
  if(!particle_config$exact & is.null(particle_config$num_mcmc)){
    num_mcmc <- ceiling(model_config$N * log(model_config$N));
    warning("Running nlogn iterations of MCMC for CondBern")
  }
  if(particle_config$verbose) cat("[ess threhold is", particle_config$ess_threshold,"]\n");
  logweights <- logw <- rep(0, particle_config$num_particles)
  logW <- rep(log(1 / particle_config$num_particles), particle_config$num_particles);
  ### logw: cumulative weights, used for resampling
  ### logW: log of normalized logw
  ### logweights: potential at the current step, use these to compute normalizing constant
  ess <- rep(0, num_observations);
  ess[1] <- 1;
  ancestors <- c(1:particle_config$num_particles)
  if(particle_config$clock){
    runtimes <- rep(NA,num_observations)
    startstopwatch <- as.numeric(Sys.time())
  }
  # treat t = 0 differently 
  current_support_y <- y[0 + 1] : model_config$N
  if (all(model_config$alpha0 == model_config$alpha0[1])){
    if(particle_config$verbose) cat("homogenous intial infection probability \n")
    logvt <- dbinom(x = current_support_y, size = model_config$N , prob = model_config$alpha0[1], log = TRUE) + dbinom(x = y[0 + 1], size = current_support_y, prob = model_config$rho, log = T)
  }else{
    if(particle_config$exact){
      logvt <- logdpoisbinom_cpp(alpha = model_config$alpha0)[current_support_y + 1] + dbinom(x = y[1], size = current_support_y, prob = model_config$rho, log = T)  
    }else{
      logvt <-  logdtranspoisson_approx(alpha = model_config$alpha0, y = current_support_y) + dbinom(x = y[1], size = current_support_y, prob = model_config$rho, log = T)  ;
    }
    logvt[is.na(logvt)] <- -Inf ## for the numerical problems in poisbinom package 
  }
  loglikelihood[0 + 1] <- lw.logsum(logvt) ## p(y0)
  if(particle_config$clock) runtimes[1] <- Sys.time()[1] 
  ## normalize to get conditional density of i0 | y0 
  vt <- lw.normalize(logvt)
  ## sample i0 | y0
  its <- sample(current_support_y, size = particle_config$num_particles, replace = T, prob = vt)
  ## sample x0 | i0
  if (all(model_config$alpha0 == model_config$alpha0[1])){
    if(particle_config$verbose) cat("homogenous intial infection probability \n")
    xts <- matrix(0, nrow = model_config$N, ncol = particle_config$num_particles)
    for (p in 1:particle_config$num_particles){
      xts[sample(x = model_config$N, size = its[p], replace = F) , p] <- 1
    }
  }else{
    xts <- sapply(its, function(ii) rcondbern(sum_x = ii, alpha = model_config$alpha0,  exact = particle_config$exact, num_sample = particle_config$num_mcmc)) ;
  }
  if(particle_config$save_particles)  particles[,1,] <- xts
  ## for t = 0, no need to resample because weights are the same
  for (t in 1:(num_observations - 1)){
    ## calculate p(x[t+1,k] = 1 | xt)
    alphats <- apply(xts, MARGIN = 2, FUN = function(xx) sis_get_alpha(agent_state = xx, model_config = model_config))
    ## calculat the weights p(y[t+1] | xt) and conditional pmf p(i[t+1]| xt, y[t+1]) 
    current_support_y <- y[t + 1] : model_config$N
    logvt <- vt <- matrix(data = NA, nrow = length(current_support_y), ncol = particle_config$num_particles)
    logweights <- rep(NA, particle_config$num_particles)
    for (p in 1 : particle_config$num_particles){
      if(particle_config$exact){
        logvt[,p] <- logdpoisbinom_cpp(alpha = alphats[,p])[current_support_y + 1] + dbinom(x = y[t + 1], size = current_support_y, prob = model_config$rho, log = TRUE)
      }else{
        logvt[,p] <- logdtranspoisson_approx(alpha = alphats[,p], y = current_support_y) +  dbinom(x = y[t + 1], size = current_support_y, prob = model_config$rho, log = TRUE);
      }
      logweights[p] <- lw.logsum(logvt[,p]) 
      vt[,p] <- lw.normalize(logvt[,p])
    }
    logweights[is.na(logweights)] <- - Inf
    loglikelihood[t + 1] <- lw.logsum(logweights) - log(particle_config$num_particles)
    if(particle_config$verbose) cat("t = " , t ," loglikelihood = ", loglikelihood[t+1], " ")
    if(all(is.infinite(logweights))){## if particles degenerate
      break
    }
    if(particle_config$clock) runtimes[t+1] <- Sys.time()
    ## resample: normalize the log weights
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    ess[t + 1] <- 1 / sum(weights**2) / particle_config$num_particles;
    ## resample steps 
    if(ess[t + 1] < particle_config$ess_threshold){
      ancestors <- sample.int(n = particle_config$num_particles, prob = weights, replace = TRUE);
      logW <- rep(-log(particle_config$num_particles), particle_config$num_particles);
      logw <- rep(0, particle_config$num_particles);
      if(particle_config$verbose) cat("ess at step t =", t + 1, "is", ess[t + 1],"\nRESAMPLE!\n");
    }else{
      ancestors <- 1 : particle_config$num_particles;
      logW <- log(weights);
      if(particle_config$verbose) cat("ess at step t =", t + 1, "is", ess[t + 1],"\n");
    }
    ## sample it given y and sample xt given it
    if (y[t + 1] == model_config$N){
      its <- rep(N, particle_config$num_particles)
      xts <- matrix(1, nrow = model_config$N, ncol = particle_config$num_particles)
    }else{
      # sample its with ancestors
      its <- sapply(1:particle_config$num_particles, function(iparticle) sample(x = current_support_y, size = 1, prob = vt[,ancestors[iparticle]]));
      xts <- sapply(1:particle_config$num_particles, function(iparticle) rcondbern(sum_x = its[iparticle], alpha = alphats[,ancestors[iparticle]], exact = particle_config$exact, num_sample = particle_config$num_mcmc));
    }
    if(particle_config$save_genealogy) eves <- eves[ancestors]
    if(particle_config$save_particles){
      particles[,,] <- particles[, ,ancestors]
      particles[,t+1,] <- xts
    }
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
