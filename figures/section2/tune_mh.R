## this script define the mh algorithm and gibbs algorithm for inference in the static model
## if tune == TRUE, it will also tune the stepsize for 

rm(list = ls())
library(agents)
set.seed(2020)

## load observations
load("~/Dropbox/AgentBasedModels/agents/figures/section2/data_static_observations.RData")
save_image <- FALSE

## run code to tune MH kernel stepsize and num_particles for IS
tune <- FALSE

## importance sampling settings
num_particles <- 20
particle_config <- list(num_particles = num_particles)

## simulation settings 
if (tune){
  # num_burn <- 5000
  # num_mcmc <- 50000 + num_burn
  num_burn <- 500
  num_mcmc <- 5000 + num_burn
  step_size <- 0.2
}

# choose the right num_sample for importance sampling
# [MSE of log(z_hat) is  0.6061536 with  20 particles]
# [with step size 0.2 ,MH acceptance rate is  0.2397818 ]
if(tune){
  exact_loglik <- static_loglikelihood_marginal(y = y[size_N], parameters = dgp, model_config = config_list[[size_N]], method = "exact")
  is_loglik <- sapply(1:1000,function(x) static_loglikelihood_marginal(y = y[size_N], parameters = dgp, model_config = config_list[[size_N]], method = "mc", particle_config = list(num_particles = num_particles)))
  cat("[MSE of log(z_hat) is ", mean((exact_loglik - is_loglik)**2), "with ", num_particles, "particles]\n")
}

## prior
lprior <- function(param){
  ## compute the log prior density
  ## param is [beta , rho]
  ## beta in R and rho in (0,1]
  ## standard normal prior for beta and uniform prior on rho 
  dnorm(param$beta, log = T)
  ## the prior on rho; uniform (0,1)
}


## proposal kernel with Jacobian adjustments
## proposal on rho is such that all values are in [0,1]
## thus we don't need to check that the proposed rho is in [0,1]
q <- function(param){
  ## returns proposed parameter 
  z <- step_size * rnorm(2)
  proposal <- param 
  proposal$beta <- param$beta + z[1]
  proposal$rho <- expit(logit(param$rho) + z[2])
  return(proposal)
}

## log Jacobian to account for proposal distribution on rho?
ljacobian <- function(param){
  + log(param$rho) + log(1 - param$rho)
}


## mh chain
mh <- function(num_mcmc, llik){
  ## store chain for beta and rho
  b_chain <- rep(NA, num_mcmc)
  r_chain <- rep(NA, num_mcmc)
  ## store acceptance indicator
  accept <- rep(0, num_mcmc)
  ## record wall clock time
  start_time <- Sys.time()
  ## initialize Markov chain and 
  ## store state of Markov chain in a list
  param <- list(beta = rnorm(1), rho = runif(1))
  ## compute posterior and Jacobian term of the proposal at current state
  lpost_current <- llik(param) + lprior(param) + ljacobian(param)
  ## if current posterior is -Infinity, set new initial point
  while(is.infinite(lpost_current)){
    param <- list(beta = rnorm(1), rho = runif(1))
    lpost_current <- llik(param) + lprior(param) + ljacobian(param)
  }
  ## at each iteration
  for (iter in 1:num_mcmc){
    ## propose new state
    proposal <- q(param)
    ## evaluate target density at proposal + log Jacobian
    lpost_proposal <- lprior(proposal) + llik(proposal) + ljacobian(proposal)
    ## MH acceptance step
    if (log(runif(1)) < (lpost_proposal - lpost_current)){
      accept[iter] <- 1
      param <- proposal
      lpost_current <- lpost_proposal
    }
    ## record current state
    b_chain[iter] <- param$beta
    r_chain[iter] <- param$rho
  }
  ## record wall clock time
  end_time <- Sys.time()
  ## return chains, run-time, acceptances
  return(list(beta_chain = b_chain, rho_chain = r_chain, 
              elapse = as.numeric(end_time - start_time, units = "secs"), 
              acceptance = accept))
}

## choose the right step-size for kernel q using the largest population
if(tune){
  llik <- function(param){
    static_loglikelihood_marginal(y = y[size_N], parameters = param, model_config = config_list[[size_N]], method = "exact")
  }
  # llik(list(beta = rnorm(1), rho = runif(1)))
  result <- mh(num_mcmc, llik)
  cat("[with step size", step_size, "MH acceptance rate is", round(100*mean(result$acceptance),1),"%]\n")
  cat("[", num_mcmc, "MH iterations performed in ", round(result$elapse,1), "seconds]\n")
  # matplot(result$beta_chain, type = 'l')
  # matplot(result$rho_chain, type = 'l')
  # plot(result$beta_chain, result$rho_chain)
  # plot(result$beta_chain, logit(result$rho_chain))
}

## Gibbs sampler, alternating between updates of the latent individuals 'x'
## and updates of the parameters given the latent individuals

gibbs <- function(num_mcmc, llik_complete, xx_gibbs, xx_swap, population_size){
  ## store parameter chains
  b_chain <- rep(NA, num_mcmc)
  r_chain <- rep(NA, num_mcmc)
  accept <- rep(0, num_mcmc)
  ## record wall-clock time
  start_time <- Sys.time()
  ## initialize parameter
  param <- list(beta = rnorm(1), rho  = runif(1))
  ## initialize latent individuals... all ones??
  xx <- rep(1, population_size)
  ## evaluate current posterior density
  lpost_current <- llik_complete(param, xx) + lprior(param) + ljacobian(param)
  while(is.infinite(lpost_current)){
    ##  if target density is -Infinity, find new starting point
    param <- list(beta = rnorm(1), rho  = runif(1))
    lpost_current <- llik_complete(param, xx) + lprior(param) + ljacobian(param)
  }
  ## for each iteration
  for (iter in 1:num_mcmc){
    # first update parameters given x  
    proposal <- q(param)
    lpost_proposal <- llik_complete(proposal, xx) + lprior(proposal) + ljacobian(proposal)
    if (log(runif(1)) < (lpost_proposal - lpost_current)){
      accept[iter] <- 1
      param <- proposal
    }
    b_chain[iter] <- param$beta
    r_chain[iter] <- param$rho
    ## then update x given parameters
    ## with probability 0.5 perform a systematic-scan gibbs update or a blocked gibbs update
    ## with probability 0.5 perform swaps
    if (runif(1) < 0.5){
      xx <- xx_gibbs(param, xx)
    }else{
      xx <- xx_swap(param, xx)
    }
    ## recompute target density at current parameter, given new x
    lpost_current <- llik_complete(param, xx) + lprior(param) + ljacobian(param) 
  }
  ## record wall-clock time
  end_time <- Sys.time()
  return(list(beta_chain = b_chain, rho_chain = r_chain, acceptance = accept,
              elapse = as.numeric(end_time - start_time, units = "secs")))
}

if (save_image){
  save.image(file = 'figures/section2/data_tune_mh.RData')
}



