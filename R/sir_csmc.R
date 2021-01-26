#' @title controlled SMC for SIR model
#' @description marginal likelihood estimations and samples from the posterior distribution of agent states
#' @param y , observations, length (T+1)
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @param particle_config a list specifying the particle filter, must pass the test of check_particle_config;
#' @return the particle filter outcome, must contain estimate of log marginal likelihood and effective sample size.
#' @export

sir_csmc <- function(y, model_config, particle_config, policy = NULL){
  if (is.null(policy)) policy <- sir_backward_information_filter_sumbin(y, model_config)
    ## policy[it + 1,rt + 1,t + 1] = logphi(it,rt | y[t:T])
  num_observations <- length(y)
  lognormalisingconstant <- rep(- Inf, num_observations)
  logweights <- rep(-Inf , particle_config$num_particles)
  ess <- rep(0, num_observations)
  eves <- c(1 : particle_config$num_particles)
  particles <- array(NA, dim = c(N, num_observations, particle_config$num_particles))
  it <- rt <- rep(NA, particle_config$num_particles)
  xts <- matrix(NA, ncol = particle_config$num_particles, nrow = N)
  twist_q <- normalize_twist_q <- array(NA, dim = c(N + 1, N + 1, particle_config$num_particles))
  logcondexp <- rep(-Inf, particle_config$num_particles)
  ## start the clock
  if(particle_config$clock){
    runtimes <- rep(NA,num_observations - 1)
    startstopwatch <- as.numeric(Sys.time());
  }
  # treat t = 0 differently, because alpha0
  ## compute the posterior of sum normalize to get the weights and resample to get the sum 
  ## mu0(i0,r0) * policy_0(i0,r0) for i0 = y0 : N 
  current_infection_support <- y[1] : model_config$N
  logposterior <- logdpoisbinom(model_config$alpha0)[current_infection_support + 1] +  policy[current_infection_support + 1, 1 + 0, 1] ## at t = 0, rt = 0 by alpha0 
  maxlogposterior <- max(logposterior)
  logposterior <- logposterior - maxlogposterior
  logz0 <- log(sum(exp(logposterior))) + maxlogposterior ## this is mu0(policy0)
  lognormalisingconstant[1] <- logz0 ## everything until this point looks correct
  # print(logz0)
  if (particle_config$clock) runtimes[1] <- as.numeric(Sys.time());
  normalizedposterior <- exp(logposterior) / sum(exp(logposterior)); ## use this to sample it at t = 0 
  for (particle in 1:particle_config$num_particles){## can vectorize this later 
    ### sample i0
    it[particle] <- sample(x = current_infection_support, size = 1, replace = FALSE, prob = normalizedposterior)
    rt[particle] <- 0 
    ### sample x0 given i0 and alpha0
    xts[,particle] <- rcondbern(sum_x = it[particle], alpha = model_config$alpha0, exact = TRUE)
  }
  ## now that we have samples from the initial distribution of the twisted path measure, we need to weight them by the twisted potential
  for (t in 1 : (num_observations - 1)){
    ## compute the log-weights
    ## first compute the sampling kernel for next step F[t+1](policy[t+1])(xt), which is the sum of product of policy(.) and kernel Q(x0,.)
    logcondexp <- rep(- Inf, particle_config$num_particles);
    for (p in 1 : particle_config$num_particles){
      twist_q[ , , p] <- sir_logdpoismulti_given_xxprev(xxprev = xts[, p], model_config = model_config) + policy[ ,  , t + 1];
      logcondexp[p] <- lw.logsum(twist_q[ , , p]);
    }
    ## potential = p(yt | xt) * E(policy |xt) / policy(xt)
    for (p in 1 : particle_config$num_particles){
      logweights[p] <- logcondexp[p] + dbinom(x = y[t], size = it[p], prob = model_config$rho, log = TRUE) - policy[it[p] + 1 , rt[p] + 1, t];
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
    logweights <- logweights - maxlogweights
    unnormalisedweights <- exp(logweights)
    weights <- unnormalisedweights / sum(unnormalisedweights)
    ess[t] <- 1 / sum(weights ** 2)
    lognormalisingconstant[t + 1] <- log(sum(unnormalisedweights)) + maxlogweights - log(particle_config$num_particles);
    if (particle_config$clock) runtimes[t + 1] <- as.numeric(Sys.time())
    ## resample
    ancester <- sample.int(n = particle_config$num_particles, size = particle_config$num_particles, replace = TRUE, prob = weights)
    eves <- eves[ancester]
    xts <- xts[ , ancester]
    it <- it[ ancester ]
    rt <- rt[ ancester ] 
    twist_q <- twist_q[ , , ancester]
    ## sample x[t+1] according to the twisted kernel 
    xtnext <- matrix(NA, nrow = model_config$N, ncol = particle_config$num_particles)
    for (p in 1 : particle_config$num_particles){
      ## sample it and rt 
      it_logweights <- apply(twist_q[ , , p], 1 , lw.logsum) 
      it[p] <- sample(0 : N, size = 1, prob = lw.normalize(it_logweights))
      rt_weights <- lw.normalize(twist_q[it[p] + 1, , p])
      rt[p] <- sample(0 : N, size = 1, prob = rt_weights)
      ## sample x
      xtnext[ , p ] <- sir_kernel_twisted_ir(xxprev = xts[,p], it = it[p], rt = rt[p], model_config = model_config)
    } 
    if (any(it < y[t + 1])) print('it problem')
    if (any(rt +  it > N)) print('rt problem')
    xts <- xtnext
  }
  # print(lognormalisingconstant)
  # print(sum(lognormalisingconstant))
  ## need to sample xT at terminal time
  ess[num_observations] <- particle_config$num_particles
  result <- list(lognormalisingconstant = lognormalisingconstant, logz = sum(lognormalisingconstant), ess = ess[-1] / particle_config$num_particles, eves = eves)
  if (particle_config$clock){
    runtimes <- runtimes - startstopwatch
    result[["runtimes"]] <- runtimes
    result[["totaltime"]] <- sum(runtimes)
  }
  return(result)
}

