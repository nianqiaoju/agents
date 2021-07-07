#' @title controlled SMC sampler for the boarding school dataset on complete networks.
#' @description marginal likelihood estimations. 
#' @param y , observations, length (T+1)
#' @param N, population size
#' @param alpha0, double, initial infection probability
#' @param lambda, double, transmission potential
#' @param gamma, double, recovery rate
#' @param rho, double, reporting rate 
#' @param logpolicy a 3-d array storing the bif information filter
#' @param num_particles number of particles, default to 20
#' @param ess_threshold adaptive resampling ESS threshold, default to 0.5
#' @return estimated log marginal likelihood
#' @export


boarding_csmc_full_fast <- function(y, N, alpha0, lambda, gamma,rho, logpolicy, num_particles = 20, ess_threshold = 0.5){
  ## logpolicy[st + 1, it + 1,t + 1] = logpsi_t(st,it)
  num_observations <- length(y);
  lognormalisingconstant <- rep(- Inf, num_observations); ## sum of lognormalising constant is final log marginal likelihood estimate
  ## logw: cumulative weights, used for resampling
  ## logW: log of normalized logw
  logweights <- logw  <- logW <- rep(0, num_particles); ## average of weights is one step normalisingconstant
  ## storage and samples
  stitindex <- it <- st <- integer(num_particles);
  xts <- matrix(0, ncol = num_particles, nrow = N);
  alphas2i <- alphai2i <-  matrix(0, nrow = N, ncol = num_particles);
  logf <- fpsi <- matrix(-Inf, nrow = dim(logpolicy)[1], ncol = num_particles); 
  ## logf stores the kernel f(snext,inext given xnow)
  ## fpsi this is the product of f * psi 
  logcondexp <- rep(-Inf, num_particles); ## this is the log normalizing constant of fpsi
  t <- 0; ## treat t = 0 differently, because r0 = 0 by alpha0
  ## mu0(s0,i0) * policy_0(s0,i0) for (s0,i0) in supp(s0,i0);
  ## the normalising constant is logz0
  ldx0 <- function(si){
    if(si[2] + si[1] == N){
      return(dbinom(x = si[2], size = N, prob = alpha0, log = TRUE))
    }else{
      return(-Inf)
    }
  }
  logposterior <- apply(lowdimstates, 1, function(si) ldx0(si));
  logposterior <- logposterior + logpolicy[, t + 1];
  maxlogposterior <- max(logposterior);
  logposterior <- logposterior - maxlogposterior;
  logz0 <- log(sum(exp(logposterior))) + maxlogposterior; ## this is mu0(policy0)
  lognormalisingconstant[t + 1] <- logz0;
  normalizedposterior <- exp(logposterior) / sum(exp(logposterior)); ## use this to sample it at t = 0 
  stitindex <- sample.int(dim(lowdimstates)[1], size = num_particles, replace = TRUE, prob = normalizedposterior);
  it <- lowdimstates[stitindex, 2];
  st <- lowdimstates[stitindex, 1];
  for (p in 1:num_particles){
    xts[sample.int(n = N, size = it[p], replace = FALSE), p] <- 1; 
  }
  for (t in 1 : (num_observations - 1)){
    ### step 1 - compute weights of each particle
    ## update alphas2i, alphai2i, and logf
    boarding_logf_update_full(logf, alphas2i, alphai2i, xts,lambda,gamma, N);
    for (p in 1 : num_particles){
      fpsi[ , p] <- logf[,p] + logpolicy[ ,  t + 1];
      logcondexp[p] <- lw.logsum(fpsi[ , p]); ## expectation of psi[t+1] given x[t]
    }
    ## unnormalized weight = p(yt | xt) * expecation(policy[t+1] |xt) / policy[t](xt)
    for (p in 1 : num_particles){
      logweights[p] <- logcondexp[p] + dbinom(x = y[t + 1], size = it[p], prob = rho, log = TRUE) - logpolicy[stitindex[p], t];
    }
    ### step 2 -- compute normalizing constant given the weights
    maxlogweights <- max(logweights);
    if (is.infinite(maxlogweights)){
      return(-Inf);
    }
    logweights <- logweights - maxlogweights;
    lognormalisingconstant[t + 1] <- log(sum(exp(logweights))) + maxlogweights - log(num_particles);
    ### step 3 - normalize the weights and (adaptively) resample particles
    ## logw, logW and ess are cumulative
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    ess <- 1 / sum(weights ** 2) / num_particles;
    # cat("t = ", t, "with ess = ", ess, "\n");
    if(ess < ess_threshold | any(is.infinite(logweights))){## adaptive resampling
      # cat(which(weights == 0), "\n");
      ancester <- sample.int(n = num_particles, prob = weights, replace = TRUE); ## resample
      # cat(ancester, "\n");
      logw <- logW <- rep(0, num_particles); ## update logw and logW
      ## re-indexing
      xts <- xts[ , ancester];
      alphas2i <- alphas2i[, ancester];
      alphai2i <- alphai2i[, ancester];
      fpsi <- fpsi[, ancester];
      rm(ancester);
    }else{
      logW <- log(weights); ## this is cumulative.
    }
    ## step 4 -- propose new particles, first sample (s,i) and then sample x
    for (p in 1 : num_particles){
      stitindex[p] <- sample.int(dim(lowdimstates)[1], size = 1, prob = lw.normalize(fpsi[,p]));
      it[p] <- lowdimstates[stitindex[p], 2];
      st[p] <- lowdimstates[stitindex[p], 1];
    } 
    xts <- boarding_sample_x_given_si(xts, alphas2i, alphai2i, st, it, N, num_particles);
  }
  return(sum(lognormalisingconstant));
}