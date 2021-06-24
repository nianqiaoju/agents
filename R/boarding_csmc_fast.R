#' @title controlled SMC sampler for the boarding school dataset.
#' @description marginal likelihood estimations. WARNING: current implementation assumes fully connected network.
#' @param y , observations, length (T+1)
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @param logpolicy a 3-d array storing the bif information filter
#' @param num_particles number of particles, default to 20
#' @param ess_threshold adaptive resampling ESS threshold, default to 0.5
#' @return estimated log marginal likelihood
#' @export


boarding_csmc_fast <- function(y, model_config, logpolicy, num_particles = 20, ess_threshold = 0.5){
  ## logpolicy[st + 1, it + 1,t + 1] = logpsi_t(st,it)
  num_observations <- length(y);
  lognormalisingconstant <- rep(- Inf, num_observations); ## sum of lognormalising constant is final log marginal likelihood estimate
  ## logw: cumulative weights, used for resampling
  ## logW: log of normalized logw
  logweights <- logw  <- logW <- rep(0, num_particles); ## average of weights is one step normalisingconstant
  ## storage and samples
  stitindex <- it <- st <- integer(num_particles);
  xts <- matrix(NA, ncol = num_particles, nrow = model_config$N);
  alphas2i <- alphai2i <-  matrix(0, nrow = model_config$N, ncol = num_particles);
  logf <- fpsi <- matrix(-Inf, nrow = dim(logpolicy)[1], ncol = num_particles); 
  ## logf stores the kernel f(snext,inext given xnow)
  ## fpsi this is the product of f * psi 
  logcondexp <- rep(-Inf, num_particles); ## this is the log normalizing constant of fpsi
  t <- 0; ## treat t = 0 differently, because r0 = 0 by alpha0
  ## mu0(s0,i0) * policy_0(s0,i0) for (s0,i0) in supp(s0,i0);
  ## the normalising constant is logz0
  logmu <- logdpoisbinom(model_config$alpha0);
  logposterior <- apply(lowdimstates, 1, function(si) (si[2] + si[1] == model_config$N) * logmu[1 + si[2]]);
  logposterior <- logposterior + logpolicy[, t + 1];
  maxlogposterior <- max(logposterior);
  logposterior <- logposterior - maxlogposterior;
  logz0 <- log(sum(exp(logposterior))) + maxlogposterior; ## this is mu0(policy0)
  lognormalisingconstant[t + 1] <- logz0;
  normalizedposterior <- exp(logposterior) / sum(exp(logposterior)); ## use this to sample it at t = 0 
  stitindex <- sample.int(dim(lowdimstates)[1], size = num_particles, replace = TRUE, prob = normalizedposterior);
  it <- lowdimstates[stitindex, 2];
  st <- lowdimstates[stitindex, 1];
  for (p in 1:num_particles){## can vectorize this later 
    xts[,p] <- as.integer(rcondbern(sum_x = it[p], alpha = model_config$alpha0, exact = TRUE));
  }
  for (t in 1 : (num_observations - 1)){
    ## expectation of psi[t+1] given x[t]
    ## update alphas2i, alphai2i, and logf
    boarding_logf_update_sparse(logf, alphas2i, alphai2i, xts, model_config$lambda, model_config$gamma, model_config$neighbors, model_config$N);
    for (p in 1 : num_particles){
      fpsi[ , p] <- logf[,p] + logpolicy[ ,  t + 1];
      logcondexp[p] <- lw.logsum(fpsi[ , p]);
    }
    ## unnormalized weight = p(yt | xt) * expecation(policy[t+1] |xt) / policy[t](xt)
    for (p in 1 : num_particles){
      logweights[p] <- logcondexp[p] + dbinom(x = y[t + 1], size = it[p], prob = model_config$rho, log = TRUE) - logpolicy[stitindex[p], t];
    }
    ## normalise the weights
    maxlogweights <- max(logweights);
    if (is.infinite(maxlogweights)){
      return(-Inf);
    }
    logweights <- logweights - maxlogweights;
    normalizedweights <- exp(logweights) / sum(exp(logweights));
    lognormalisingconstant[t + 1] <- log(sum(exp(logweights))) + maxlogweights - log(num_particles);
    ## logw, logW and ess are cumulative
    logw <- logW + logweights;
    weights <- lw.normalize(logw);
    ess <- 1 / sum(weights ** 2) / num_particles;
    cat("t = ", t, "with ess = ", ess, "\n");
    if(ess < ess_threshold | any(is.infinite(logweights))){## adaptive resampling
      cat(which(weights == 0), "\n");
      ancester <- sample.int(n = num_particles, prob = weights, replace = TRUE); ## resample
      cat(ancester, "\n");
      logw <- logW <- rep(0, num_particles); ## update logw and logW
      ## re-indexing
      xts <- xts[ , ancester];
      alphas2i <- alphas2i[, ancester];
      alphai2i <- alphai2i[, ancester];
      fpsi <- fpsi[, ancester];
      rm(ancester);
    }else{
      logW <- log(weights); ## this is cumulative. it is equivalent to logW <- normalized(logw)
    }
    ## sample xnext according to the twisted kernel 
    for (p in 1 : num_particles){
      stitindex[p] <- sample.int(dim(lowdimstates)[1], size = 1, replace = FALSE, prob = lw.normalize(fpsi[,p]));
      it[p] <- lowdimstates[stitindex[p], 2];
      st[p] <- lowdimstates[stitindex[p], 1];
    } 
    ## sammple xt given st and it
    boarding_sample_x_given_si_sparse(xts, alphas2i, alphai2i, model_config$lambda, model_config$gamma, st, it, model_config$N, num_particles);
  }
  return(sum(lognormalisingconstant));
}