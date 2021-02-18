#' @title controlled SMC for SIR model
#' @description marginal likelihood estimations and samples from the posterior distribution of agent states. This function is the development version of sir_csmc. 
#' Its policy will be a function of (st,it).
#' @param y , observations, length (T+1)
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @param particle_config a list specifying the particle filter, must pass the test of check_particle_config;
#' @return the particle filter outcome, must contain estimate of log marginal likelihood and effective sample size.
#' @export

sir_csmc_si <- function(y, model_config, particle_config, logpolicy = NULL){
  if(is.null(particle_config$ess_threshold)) stop("Please specify ess_threshold for resampling.");
  if(particle_config$verbose) cat("[ess threhold is", particle_config$ess_threshold,"]\n");
  if(is.null(logpolicy)) logpolicy <- sir_backward_information_filter_si(y, model_config);
  ## logpolicy[st + 1, it + 1,t + 1] = logpsi(st, it | y[t:T])
  num_observations <- length(y);
  ## sum of lognormalising constant is final log marginal likelihood estimate
  lognormalisingconstant <- rep(- Inf, num_observations);
  ## average weight is one step normalisingconstant
  logweights <- logw <- rep(-Inf, particle_config$num_particles)
  logW <- rep(log(1 / particle_config$num_particles), particle_config$num_particles);
  ### logw: cumulative weights, used for resampling
  ### logW: log of normalized logw
  ### logweights: potential at the current step, use these to compute normalizing constant
  ## diagonsis and genealogy
  ess <- rep(0, num_observations)
  eves <- c(1 : particle_config$num_particles)
  ## storage and samples
  particles <- array(NA, dim = c(N, num_observations, particle_config$num_particles));
  it <- st <- rep(NA, particle_config$num_particles);
  xts <- matrix(NA, ncol = particle_config$num_particles, nrow = N);
  fpsi <- array(NA, dim = c(N + 1, N + 1, particle_config$num_particles)); ## this is the product of f * psi 
  logcondexp <- rep(-Inf, particle_config$num_particles); ## this is the log normalizing constant of fpsi
  ## start the clock
  if(particle_config$clock){
    runtimes <- rep(NA,num_observations - 1)
    startstopwatch <- as.numeric(Sys.time());
  }
  ## treat t = 0 differently
  ## mu0(i0,r0) * policy_0(i0,r0) for i0 = y0 : N 
  ## the normalising constant is logz0
  if(y[1] == model_config$N){
    logz0 <- logpolicy[0 + 1, model_config$N + 1, 1] + sum(log(model_config$alpha0));
    lognormalisingconstant[1] <- logz0
    xts <- matrix(1, ncol = particle_config$num_particles, nrow = N);
    st <- rep(0, particle_config$num_particles);
    it <- rep(model_config$N, particle_config$num_particles);
  }else{
    current_infection_support <- y[1] : model_config$N;
    logposterior <- sapply(current_infection_support, function(icount) logpolicy[model_config$N - icount + 1, icount + 1, 1]) + logdpoisbinom(model_config$alpha0)[current_infection_support + 1] ## at t = 0, rt = 0 by alpha0.
    maxlogposterior <- max(logposterior);
    logposterior <- logposterior - maxlogposterior
    logz0 <- log(sum(exp(logposterior))) + maxlogposterior ## this is mu0(policy0)
    lognormalisingconstant[1] <- logz0 
    if (particle_config$clock) runtimes[1] <- as.numeric(Sys.time());
    normalizedposterior <- exp(logposterior) / sum(exp(logposterior)); ## use this to sample it at t = 0 
    for (p in 1:particle_config$num_particles){## can vectorize this later 
      ### sample i0 and s0
      it[p] <- sample(x = current_infection_support, size = 1, replace = FALSE, prob = normalizedposterior)
      st[p] <- model_config$N - it[p];
      ### sample x0 given i0 and alpha0
      xts[,p] <- rcondbern(sum_x = it[p], alpha = model_config$alpha0, exact = TRUE)
    }
  }
  ## after the initial step
  for (t in 1 : (num_observations - 1)){
    ## now that we have samples from the initial distribution of the twisted path measure, we need to weight them by the twisted potential
    ## first compute the sampling kernel for next step F[t+1](policy[t+1])(xt), which is the sum of product of psi[t]() and kernel f(x[t-1],.)
    logcondexp <- rep(- Inf, particle_config$num_particles);
    for (p in 1 : particle_config$num_particles){
      fpsi[ , , p] <- sir_logdpoismulti_given_xxprev_si(xxprev = xts[, p], model_config = model_config) + logpolicy[ ,  , t + 1];
      logcondexp[p] <- lw.logsum(fpsi[ , , p]);
    }
    ## unnormalized weight = p(yt | xt) * E(policy |xt) / policy(xt)
    for (p in 1 : particle_config$num_particles){
      logweights[p] <- logcondexp[p] + dbinom(x = y[t], size = it[p], prob = model_config$rho, log = TRUE) - logpolicy[st[p] + 1, it[p] + 1 , t];
    }
    ## normalise the weights
    maxlogweights <- max(logweights);
    if (is.infinite(maxlogweights)){
      print("break")
      result <- list(lognormalisingconstant = lognormalisingconstant, logz = sum(lognormalisingconstant), ess = ess)
      if (particle_config$clock){
        runtimes <- runtimes - startstopwatch
        result[["runtimes"]] <- runtimes
        result[["totaltime"]] <- sum(runtimes)
      }
      return(result)
    }
    logweights <- logweights - maxlogweights;
    lognormalisingconstant[t + 1] <- log(sum(exp(logweights))) + maxlogweights - log(particle_config$num_particles);
    if (particle_config$clock) runtimes[t + 1] <- as.numeric(Sys.time())
    ## logw, logW and ess are cumulative
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    ess[t] <- 1 / sum(weights ** 2)
    ## adaptive resampling
    if(ess[t] < particle_config$ess_threshold){
      ancester <- sample.int(n = particle_config$num_particles, prob = weights, replace = TRUE);
      logW <- rep(-log(particle_config$num_particles), particle_config$num_particles);
      logw <- rep(0, particle_config$num_particles);
      if(particle_config$verbose) cat("RESAMPLE!\ness at step t =", t, "is", ess[t],"\n");
    }else{
      ancester <- 1 : particle_config$num_particles;
      logW <- log(weights);
      if(particle_config$verbose) cat("ess at step t =", t, "is", ess[t],"\n");
    }
    eves <- eves[ancester]
    xts <- xts[ , ancester]
    it <- it[ ancester ]
    st <- st[ ancester ] 
    fpsi <- fpsi[ , , ancester]
    ## sample x[t+1] according to the twisted kernel 
    xtnext <- matrix(NA, nrow = model_config$N, ncol = particle_config$num_particles)
    for (p in 1 : particle_config$num_particles){
      ## sample it and st 
      it_logweights <- apply(fpsi[ , , p], 2 , lw.logsum) ## the column indices correspond to icount
      it[p] <- sample(0 : N, size = 1, prob = lw.normalize(it_logweights))
      st_weights <- lw.normalize(fpsi[, it[p] + 1, p])
      st[p] <- sample(0 : N, size = 1, prob = st_weights)
      ## sample x given s and i 
      xtnext[ , p ] <- sir_kernel_twisted_si(xxprev = xts[,p], st = st[p], it = it[p], model_config = model_config)
    } 
    if (any(it < y[t + 1])) print('it problem')
    if (any(st +  it > N)) print('st problem')
    xts <- xtnext
  }
  ## terminal time t = T, remember one sample from xT
  t <- num_observations; ## all the weights are 1
  xT <- xts[, 1];
  # print(lognormalisingconstant)
  # print(sum(lognormalisingconstant))
  ess[num_observations] <- particle_config$num_particles;
  result <- list(lognormalisingconstant = lognormalisingconstant, 
                 logz = sum(lognormalisingconstant), 
                 xT = xT,
                 ess = ess[-1] / particle_config$num_particles, 
                 eves = eves)
  if (particle_config$clock){
    runtimes <- runtimes - startstopwatch;
    result[["runtimes"]] <- runtimes;
    result[["totaltime"]] <- sum(runtimes);
  }
  return(result);
}

