### given the dataset generated in setup_hetero
### we would like to perform parameter inference


### this script contains helper functions to 
### (1) update model_config given parameters 
### (2) random walk kernel on parameters
### (3) defines prior distribution
### (4) pmcmc algorithm to impute parameters

### these functions are specific to d = 2 (each agent has only 1 feature)

### update model_config given parameters
### each model_config must contain the following items:
### N, adjacency_matrix_b, network_type, alpha0, lambda, gamma, rho
### if model_config is intended for cSMC, then we also need policy = sis_backward_information_filter_sumbin(y, policy)

update_model_config <- function(model_config, parameters){
  model_config$alpha0 <- get_lambda(parameters[1], parameters[2]);
  model_config$lambda <- get_lambda(parameters[3], parameters[4]);
  model_config$gamma <- get_lambda(parameters[5],parameters[6]);
  model_config$rho <- parameters[7];
  if(! is.null(model_config$policy)) model_config$policy <- sis_backward_information_filter_sumbin(y, model_config);
  return(model_config);
}


### random walk kernel on parameters
### need to transfer rho to logit scale
qkernel <- function(parameters, step_size){
  z <- step_size * rnorm(7);
  proposal <- parameters;
  proposal[1:6] <- parameters[1:6] + z[1:6];
  proposal[7] <- expit(logit(parameters[7]) + z[7]);
  return(proposal)
}

ljacobian <- function(parameters){
  + log(parameters[7]) + log(1 - parameters[7]);
}

### prior
lprior <- function(parameters){
  ## compute the log prior density
  ## iid standard normal prior for each beta parameter
  ## uniform prior on rho 
  sum(dnorm(parameters[1:6], sd = 3, log = T))
}

csmc <- function(model_config){
  sis_csmc(y, model_config, particle_config)$log_final_likelihood;
}

bpf <- function(model_config){
  sis_bpf(y, model_config, particle_config)$log_final_likelihood;
}

pmcmc <- function(num_mcmc, pf, init_config, step_size){
  ### pf is a function that takes model_config as an input and returns logz 
  ### save and output these values
  chain <- matrix(nrow = num_mcmc, ncol = 7); ## each row is a parameter
  accept_chain <- rep(0, num_mcmc);
  ### instantiate parameters
  curr_param <- c(rnorm(6), runif(1));
  curr_config <- update_model_config(init_config, curr_param);
  curr_lpost <- pf(curr_config) + lprior(curr_param) + ljacobian(curr_param);
  if(is.infinite(curr_lpost)) stop("initialization of the parameters get 0 likelihood!");
  # prop_config <- curr_config;
  for(imcmc in 1 : num_mcmc){
    ## propose
    prop_param <- qkernel(curr_param, step_size);
    prop_config <- update_model_config(curr_config, prop_param);
    prop_lpost <- pf(prop_config) + lprior(prop_param) + ljacobian(prop_param);
    ## accept and reject
    if (log(runif(1)) < (prop_lpost - curr_lpost)){
      accept_chain[imcmc] <- 1
      curr_param <- prop_param;
      curr_lpost <- prop_lpost;
    }
    ## save the parameters
    chain[imcmc, ] <- curr_param;
    if(imcmc %% 100 ==0) cat("[",imcmc, "out of", num_mcmc, "iterations]\n");
  }
  return(list(parameters_chain = chain, accept_chain = accept_chain));
}


## TO DO: add gibbs sampler wrappers here

tune_ <- FALSE
if(tune_){
  cat("WARNING: this is running code to tune step size!")
  load("figures/section3/data_setup_hetero.RData");
  dgp_parameters <- c(log(1 - 1/N), 0, dgp_bl1, dgp_bl2, dgp_bg1, dgp_bg2, dgp_config$rho);
  library(agents)
  particle_config <- list(ess_threshold = 0.5,
                          num_particles = 64,
                          clock = FALSE,
                          save_genealogy = FALSE,
                          save_particles = FALSE,
                          verbose = FALSE,
                          exact = TRUE)
  pf <- csmc;
  dgp_config$policy <- sis_backward_information_filter_sumbin(y, dgp_config);
  num_mcmc <- 10000 ;
  init_config <- dgp_config;
  step_size <- 0.2;
  set.seed(2020);
  result <- pmcmc(num_mcmc, pf, init_config, step_size);
  cat(mean(result$accept_chain), "\n");
  save.image("figures/section3/data_tune_pmcmc.R");
  
  plot.new()
  matplot(result$parameters_chain, col = 1: 7, type = "l", lty = 1)
  legend("topleft", paste(c(1:7)), col= c(1:7) , cex= 0.8,fill= c(1:7))
  abline(h = dgp_parameters, col =  1:7);
  
  colMeans(result$parameters_chain[5000:10000, ])
  print(dgp_parameters)
}