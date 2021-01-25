#' @title controlled SMC sampler for SIS model with population-level observations
#' @description controlled SMC sampler for SIS model with population-level observations
#' @param y a vector of length (T+1)
#' @param model_config a list containing:
#' @param particle_config 
#' @return A list containing 
#' \itemize{
#' \item 'log_final_likelihood' : log marginal likelihood estimator;
#' \item 'xT': sample from xT given y(0:T);
#' \item 'particles': [N, (T+1), num_particles] array storing the particles, if particle_config$save_particles = TRUE'; 
#' \item 'eves' : the unique ancesters of each particle in the genealogy structure, if save_genealogy = TRUE;
#' \item 'ess' : effective sample size at each step
#' \item 'runtime' : elapse time at each step, if clock = TRUE;
#' \item 'totaltime' : elapse time of running the whole particle filter, if clock = TRUE;
#' }
#' @export

sis_csmc <- function(y, model_config, particle_config){
  num_observations <- length(y)
  lognormalisingconstant <- rep(-Inf, num_observations - 1)
  logweights <- logw <- rep(0, particle_config$num_particles)
  logW <- rep(log(1 / particle_config$num_particles), particle_config$num_particles);
  ### logw: cumulative weights, used for resampling
  ### logW: log of normalized logw
  ### logweights: potential at the current step, use these to compute normalizing constant
  ess <- rep(1, num_observations);
  ancestors <- c(1:particle_config$num_particles)
  if(is.null(model_config$policy)){
    warning("policy is not given and cSMC is computing the backward information filter.\n");
    model_config$policy <- sis_backward_information_filter_sumbin(y, model_config);
  }
  if(is.null(particle_config$exact)) stop("Please specify methods to use for density evalutions and for sampling.")
  if(!particle_config$exact & is.null(particle_config$num_mcmc)){
    particle_config$num_mcmc <- ceiling(model_config$N * log(model_config$N));
    warning("Running nlogn iterations of MCMC for CondBern")
  }
  if(particle_config$save_particles)  particles <- array(NA, dim = c(model_config$N, num_observations, particle_config$num_particles))
  if(particle_config$save_genealogy){
    eves <- c(1 : particle_config$num_particles)
  } 
  if(particle_config$clock){
    runtimes <- rep(NA,num_observations - 1)
    startstopwatch <- as.numeric(Sys.time())
  }
  if(is.null(particle_config$ess_threshold)) stop("Please specify ess_threshold for resampling.");
  if(particle_config$verbose) cat("[ess threhold is", particle_config$ess_threshold,"]\n");
  ## treat t = 0 differently 
  current_support_y <- y[0 + 1] : model_config$N
  if (all(model_config$alpha0 == model_config$alpha0[1])){
    if(particle_config$verbose) cat("homogenous intial infection probability \n")
    logvt <- dbinom(x = current_support_y, size = model_config$N , prob = model_config$alpha0[1], log = TRUE) + model_config$policy[current_support_y + 1,0 + 1]
  }else{
    if(particle_config$exact){
      logvt <- logdpoisbinom_cpp(alpha = model_config$alpha0)[current_support_y + 1] + model_config$policy[current_support_y + 1,0 + 1]
    }else{
      logvt <- logdtranspoisson_approx(alpha = model_config$alpha0, y = current_support_y) + model_config$policy[current_support_y + 1,0 + 1];
    }
  }
  if(y[1] == model_config$N){
    its <- rep(y[1], particle_config$num_particles);
    xts <-  sapply(its, function(it) rcondbern(sum_x = it, alpha = model_config$alpha0, exact = particle_config$exact, num_sample = particle_config$num_mcmc))
  }else{
    ## sample i0 for each particle
    its <- sample(x = current_support_y, size = particle_config$num_particles, replace = TRUE, prob = lw.normalize(logvt))
    ## sample x0 for each particle
    xts <- sapply(its, function(it) rcondbern(sum_x = it, alpha = model_config$alpha0, exact = particle_config$exact, num_sample = particle_config$num_mcmc))
  }

  ## compute alpha for the samples
  alphats <- apply(xts, 2, function(xx) sis_get_alpha(agent_state = xx, model_config = model_config))
  ## compute weights
  first_cond_expecation <- lw.logsum(logvt)
  ## computing the conditional expectations for the next step
  current_support_y <- y[0 + 1 + 1] : model_config$N
  logvt <- matrix(NA, nrow = length(current_support_y), ncol = particle_config$num_particles)
  for(iparticle in 1 : particle_config$num_particles){
    if(particle_config$exact){
      logvt[,iparticle] <- logdpoisbinom_cpp(alpha = alphats[, iparticle])[current_support_y + 1] + model_config$policy[current_support_y  + 1, 0 + 1 + 1]
    }else{
      logvt[,iparticle] <- logdtranspoisson_approx(alpha = alphats[, iparticle], y = current_support_y) + model_config$policy[current_support_y  + 1, 0 + 1 + 1]
    }
  }
  logweights <- rep(NA, particle_config$num_particles)
  for(iparticle in 1 : particle_config$num_particles){
    logweights[iparticle] <- first_cond_expecation + 
      dbinom(x = y[0 + 1], size = its[iparticle], prob = model_config$rho, log = TRUE) +
      lw.logsum(logvt[,iparticle]) - 
      model_config$policy[its[iparticle] + 1, 0 + 1]
  }
  for (t in 1 : (num_observations - 2)){
    ## store the lognormalising constants
    lognormalisingconstant[t] <- lw.logsum(logweights) - log(particle_config$num_particles)
    ## update weights 
    logw <- logweights + logW;
    ## normalize weights
    weights <- lw.normalize(logw);
    ## effective sample size
    ess[t] <- (1 / sum(weights**2)) / particle_config$num_particles;
    if(particle_config$clock) runtimes[t] <- Sys.time()
    ## resample steps 
    if((ess[t] < particle_config$ess_threshold)){
      ancestors <- sample.int(n = particle_config$num_particles, prob = weights, replace = TRUE)
      logW <- rep(-log(particle_config$num_particles), particle_config$num_particles);
      logw <- rep(0, particle_config$num_particles);
      if(particle_config$verbose) cat("RESAMPLE!\ness at step t =", t, "is", ess[t],"\n");
    }else{
      ancestors <- 1 : particle_config$num_particles;
      logW <- log(weights);
      if(particle_config$verbose) cat("ess at step t =", t, "is", ess[t],"\n");
    }
    if(particle_config$save_genealogy) eves <- eves[ancestors]
    ## sample it and xt
    if(y[t+1] == model_config$N){
      its <- rep(model_config$N, particle_config$num_particles)
      xts <- matrix(1, nrow = model_config$N, ncol = particle_config$num_particles)
    }else{
      its <- rep(NA, particle_config$num_particles)
      xts <- matrix(NA, nrow = model_config$N, ncol = particle_config$num_particles)
      for(iparticle in 1: particle_config$num_particles){
        its[iparticle] <- sample(x = current_support_y, size = 1, replace = FALSE, prob = lw.normalize(logvt[,ancestors[iparticle]]))
        xts[,iparticle] <- rcondbern(sum_x = its[iparticle], alpha = alphats[,ancestors[iparticle]], exact = particle_config$exact, num_sample = particle_config$num_mcmc);
      }
    }
    # if(any(its != colSums(xts))){
    #   cat("its and xts does not agree")
    #   break
    # }
    ## compute alpha of the new particles
    alphats <- apply(xts, 2, function(xx) sis_get_alpha(agent_state = xx, model_config = model_config))
    ## compute new logvt to get conditional expecatations for the next step
    current_support_y <- y[t + 1 + 1] : N
    logvt <-  matrix(NA, nrow = length(current_support_y), ncol = particle_config$num_particles)
    for(iparticle in 1 : particle_config$num_particles){
      if(particle_config$exact){
        logvt[,iparticle] <- logdpoisbinom_cpp(alpha = alphats[, iparticle])[current_support_y + 1] + model_config$policy[current_support_y  + 1, t + 1 + 1];
      }else{
        logvt[,iparticle] <- logdtranspoisson_approx(alpha = alphats[, iparticle], y = current_support_y) + model_config$policy[current_support_y  + 1, t + 1 + 1];
      }
    }
    ## compute logweights
    logweights <- rep(NA, particle_config$num_particles)
    for(iparticle in 1 : particle_config$num_particles){
      logweights[iparticle] <- dbinom(x = y[t + 1], size = its[iparticle], prob = model_config$rho, log = TRUE) + 
        lw.logsum(logvt[,iparticle]) - 
        model_config$policy[its[iparticle] + 1, t + 1]
    }
    
    if(any(is.infinite(logweights))){
      cat("there is one 0 weight");
      # break
    }
    if(all(is.infinite(logweights))){
      warning("cSMC collapsed");
      result <- list(log_incremental_likelihood = lognormalisingconstant, 
                     log_final_likelihood = -Inf,
                     ess = ess);
      if(particle_config$clock){
        result$totaltime = NA;
      }
      if(particle_config$save_particles) result$particles = particles
      if(particle_config$save_genealogy) result$eves = eves
      return(result)
    }
  }
  ## treat terminal time differently
  lognormalisingconstant[num_observations - 1] <- lw.logsum(logweights) - log(particle_config$num_particles)
  logw <- logW + logweights;
  ess[num_observations - 1] <- 1 / sum(lw.normalize(logw)**2) / particle_config$num_particles;
  ## sample xT given observations : resample xts given logweights
  if(y[num_observations] == model_config$N){
    xT <- rep(TRUE, model_config$N);
  }else{
    iparticle <- sample.int(particle_config$num_particles, 1, F, lw.normalize(logweights + logW));
    iT <- sample(x = current_support_y, size = 1, replace = FALSE, prob = lw.normalize(logvt[,iparticle]))
    xT <- rcondbern(sum_x = iT, alpha = alphats[,iparticle], exact = particle_config$exact, num_sample = particle_config$num_mcmc);
  }
  if(particle_config$clock) runtimes[num_observations - 1] <- Sys.time()
  result <- list(log_incremental_likelihood = lognormalisingconstant,
                 log_final_likelihood = sum(lognormalisingconstant),
                 xT = xT,
                 ess = ess );
  if(particle_config$clock){
    runtimes <- as.numeric(runtimes - startstopwatch, units = "seconds");
    result$runtimes = runtimes
    result$totaltime = mylast_nona(runtimes);
  }
  if(particle_config$save_particles) result$particles = particles
  if(particle_config$save_genealogy) result$eves = eves
  return(result)
}
