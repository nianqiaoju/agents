#' @title controlled SMC for the smallpox dataset.
#' @description marginal likelihood estimations and samples from the posterior distribution of agent states. 
#' This function is the same as smallpox_csmc, except that it only returns the log marginal likelihood.
#' @param y , observations, length (T+1)
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @param logpolicy a 3-d array storing the bif information filter
#' @param num_particles number of particles, default to 20
#' @param ess_threshold adaptive resampling ESS threshold, default to 0.5
#' @return estimated log marginal likelihood
#' @export

smallpox_csmc_fast <- function(y, model_config, logpolicy, num_particles = 20, ess_threshold = 0.5){
  ## logpolicy[st + 1, it + 1,t + 1] = logpsi_t(st,it)
  num_observations <- length(y);
  lognormalisingconstant <- rep(- Inf, num_observations); ## sum of lognormalising constant is final log marginal likelihood estimate
  ## logw: cumulative weights, used for resampling
  ## logW: log of normalized logw
  logweights <- logw <- rep(-Inf, num_particles); ## average of weights is one step normalisingconstant
  logW <- rep(log(1 / num_particles), num_particles);
  ## storage and samples
  it <- st <- rep(NA, num_particles);
  xts <- matrix(NA, ncol = num_particles, nrow = N);
  logf <- matrix(0, nrow = model_config$N + 1, ncol = model_config$N + 1); ##stores the kernel f(snext,inext given xnow), updated for each time and each particle
  fpsi <- array(NA, dim = c(N + 1, N + 1, num_particles)); ## this is the product of f * psi 
  logcondexp <- rep(-Inf, num_particles); ## this is the log normalizing constant of fpsi
  ## treat t = 0 differently, because r0 = 0 by alpha0
  ## mu0(s0,s0) * policy_0(s0,i0) for i0 = 0 : N;
  ## the normalising constant is logz0
  current_infection_support <- 0 : model_config$N;
  logposterior <- logdpoisbinom(model_config$alpha0);
  logposterior <- logposterior + sapply(current_infection_support,function(icount) logpolicy[model_config$N - icount + 1, icount + 1, 1]);
  maxlogposterior <- max(logposterior);
  logposterior <- logposterior - maxlogposterior;
  logz0 <- log(sum(exp(logposterior))) + maxlogposterior; ## this is mu0(policy0)
  lognormalisingconstant[1] <- logz0;
  normalizedposterior <- exp(logposterior) / sum(exp(logposterior)); ## use this to sample it at t = 0 
  for (p in 1:num_particles){## can vectorize this later 
    it[p] <- sample(x = current_infection_support, size = 1, replace = FALSE, prob = normalizedposterior)
    st[p] <- model_config$N - it[p];
    xts[,p] <- rcondbern(sum_x = it[p], alpha = model_config$alpha0, exact = TRUE)
  }
  for (t in 1 : (num_observations - 1)){
    ## expectation of psi[t+1] given x[t]
    for (p in 1 : num_particles){
      sir_csmc_update_f_matrix(logf, xts[,p], model_config$lambda, model_config$gamma, model_config$N); ## changes the matrix logf
      fpsi[ , , p] <- logf + logpolicy[ ,  , t + 1];
      logcondexp[p] <- lw.logsum(fpsi[ , , p]);
    }
    ## unnormalized weight = p(yt | xt) * expecation(policy[t+1] |xt) / policy[t](xt)
    for (p in 1 : num_particles){
      logweights[p] <- logcondexp[p] + dbinom(x = y[t], size = model_config$N - it[p] - st[p], prob = model_config$rho, log = TRUE) - logpolicy[st[p] + 1, it[p] + 1 , t];
    }
    ## normalise the weights
    maxlogweights <- max(logweights);
    if (is.infinite(maxlogweights)){
      return(-Inf);
    }
    logweights <- logweights - maxlogweights;
    lognormalisingconstant[t + 1] <- log(sum(exp(logweights))) + maxlogweights - log(num_particles);
    ## logw, logW and ess are cumulative
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    ess <- 1 / sum(weights ** 2) / num_particles;
    if(ess < ess_threshold | any(is.infinite(logweights))){## adaptive resampling
      ancester <- sample.int(n = num_particles, prob = weights, replace = TRUE); ## resample
      logW <- rep(-log(num_particles), num_particles); ## update logW
      logw <- rep(0, num_particles); ## update logw
      ## re-indexing
      xts <- xts[ , ancester]
      it <- it[ ancester ]
      st <- st[ ancester ] 
      fpsi <- fpsi[ , , ancester]
    }else{
      logW <- log(weights); ## this is cumulative. it is equivalent to logW <- normalized(logw)
    }
    ## sample xnext according to the twisted kernel 
    for (p in 1 : num_particles){
      ## sample it and st 
      it_logweights <- apply(fpsi[ , , p], 2 , lw.logsum) ## the column indices correspond to icount
      it[p] <- sample(x = current_infection_support, size = 1, replace = FALSE, prob = lw.normalize(it_logweights));
      st_weights <- lw.normalize(fpsi[, it[p] + 1, p]);
      st[p] <- sample(x = current_infection_support, size = 1, replace = FALSE, prob = st_weights);
    } 
    ## sammple xt given st and it, for fully connected networks
    xts <- sir_sample_x_given_si(xts, model_config$lambda, model_config$gamma, st, it, model_config$N, num_particles);
  }
  return(sum(lognormalisingconstant));
}

