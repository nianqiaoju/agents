sir_apf <- function(y, model_config, particle_config){
	num_observations <- length(y)
  loglikelihood <- rep(-Inf,num_observations)
  if(particle_config$save_particles)  particles <- array(NA, dim = c(model_config$N, num_observations, particle_config$num_particles))
  if(particle_config$save_genealogy){
    eves <- c(1 : particle_config$num_particles)
  } 
  ess <- rep(particle_config$num_particles, num_observations)
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
    logvt <- logdpoisbinom_cpp(alpha = model_config$alpha0)[current_support_y + 1] + dbinom(x = y[1], size = current_support_y, prob = model_config$rho, log = T)  
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
    xts <- sapply(its, function(ii) rcondbern(sum_x = ii, alpha = model_config$alpha0)) 
  }
  if(particle_config$save_particles)  particles[,1,] <- xts
  ## for t = 0, no need to resample because weights are the same
  for (t in 1:(num_observations - 1)){
    ## calculate p(x[t+1,k] = 1 | xt)
    alphats <- apply(xts, MARGIN = 2, FUN = function(xx) sir_get_alpha_i(xx,model_config))
    ## calculat the weights p(y[t+1] | xt) and conditional pmf p(i[t+1]| xt, y[t+1]) 
    current_support_y <- y[t + 1] : model_config$N
    logvt <- vt <- matrix(data = NA, nrow = length(current_support_y), ncol = particle_config$num_particles)
    logweights <- rep(NA, particle_config$num_particles)
    for (p in 1 : particle_config$num_particles){
      logvt[,p] <- logdpoisbinom_cpp(alpha = alphats[,p])[current_support_y + 1] + dbinom(x = y[t + 1], size = current_support_y, prob = model_config$rho, log = TRUE)
      logweights[p] <- lw.logsum(logvt[,p]) 
      vt[,p] <- lw.normalize(logvt[,p])
    }
    logweights[is.na(logweights)] <- - Inf
    loglikelihood[t + 1] <- lw.logsum(logweights) - log(particle_config$num_particles)
    if(all(is.infinite(logweights))){## if particles degenerate
      break
    }
    if(particle_config$clock) runtimes[t+1] <- Sys.time()
    ## resample: normalize the log weights
    weights <- lw.normalize(logweights)
    ess[t] <- 1 / sum(weights**2)
    ancestors <- sample.int(n = particle_config$num_particles, replace = TRUE, prob = weights)
    ## sample it given y and sample xt given it
    xnew <- matrix(NA, nrow = model_config$N, ncol = particle_config$num_particles)
    if (y[t + 1] == model_config$N){
      its <- rep(N, particle_config$num_particles)
      xnew <- matrix(1, nrow = model_config$N, ncol = particle_config$num_particles)
    }else{
      # sample its with ancestors
      its <- sapply(1:particle_config$num_particles, function(iparticle) sample(x = current_support_y, size = 1, prob = vt[,ancestors[iparticle]]))
      xnew <- sapply(1:particle_config$num_particles, function(iparticle)
      	sir_kernel_twisted_i(xts[,ancestors[iparticle]], its[iparticle], model_config))
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
                 ess = ess / particle_config$num_particles)
  if(particle_config$clock){
    runtimes <- runtimes - startstopwatch
    result$runtimes = runtimes
    result$totaltime = runtimes[length(runtimes)]
  }
  if(particle_config$save_particles) result$particles = particles
  if(particle_config$save_genealogy) result$eves = eves
  return(result)
}
