rm(list = ls())
library(agents)
load("figures/section3/data_setup_hetero.RData");
source("figures/section3/setup_mcmc.R");

## we will predict y((obs_length  + 1 ) : 90) given y(0:obs_length)

obs_length <- 30;
yobs <- y[1 : (obs_length  + 1 )]; ## observing y_{0:t};
## BIF based policy evaluation is based on the observed sequence y[0:t];
dgp_config$policy <- sis_backward_information_filter_sumbin(yobs , dgp_config);

## define the MH chain;
burn_in <- 5000;
num_mcmc <-  10**5 + burn_in;
dgp_parameters <- c(log(1/(N-1)), 0, dgp_bl1, dgp_bl2, dgp_bg1, dgp_bg2, dgp_config$rho);
step_size <- 0.5; ## tuned so that acceptance probability ~ 0.3
init_config <- dgp_config;


## define the particle filter used in the PMCMC
particle_config <- list(ess_threshold = 0.5,
                        num_particles = 128,
                        clock = FALSE,
                        save_genealogy = FALSE,
                        save_particles = FALSE,
                        verbose = FALSE,
                        exact = TRUE)

## redefine update_model_config
rm(update_model_config);
update_model_config <- function(model_config, parameters){
  model_config$alpha0 <- get_lambda(parameters[1], parameters[2]);
  model_config$lambda <- get_lambda(parameters[3], parameters[4]);
  model_config$gamma <- get_lambda(parameters[5],parameters[6]);
  model_config$rho <- parameters[7];
  if(! is.null(model_config$policy)) model_config$policy <- sis_backward_information_filter_sumbin(yobs, model_config); 
  ## BIF based policy evaluation is based on the observed sequence y[0:t];
  return(model_config);
}


### save and output these values
xT_chain <- matrix(nrow = num_mcmc, ncol = N); ## each row is a sample of xT
param_chain <- matrix(nrow = num_mcmc, ncol = 7); ## each row is a parameter
accept_chain <- rep(0, num_mcmc);
lpost_chain <- rep(0, num_mcmc);
### instantiate parameters
set.seed(2020);
# curr_param <- c(rnorm(6,  sd = 3), runif(1)); ### RANDOM INIT
# curr_config <- update_model_config(init_config, curr_param); ### RANDOM INIT
curr_param <- dgp_parameters; ### THIS INITIALIZES FROM DGP
curr_config <- dgp_config;
init_param <- curr_param;
print(init_param);
pf_result <- sis_csmc(yobs, curr_config, particle_config); ## this step includes a BIF based policy evaluation;
curr_lpost <- pf_result$log_final_likelihood + lprior(curr_param) + ljacobian(curr_param);
curr_xT <- pf_result$xT;
if(is.infinite(curr_lpost)) stop("initialization of the parameters get 0 likelihood!");
for(imcmc in  1 : num_mcmc){
  ## propose
  prop_param <- qkernel(curr_param, step_size);
  prop_config <- update_model_config(curr_config, prop_param); ## this step includes a BIF based policy evaluation;
  pf_result <- sis_csmc(yobs, prop_config, particle_config);
  prop_lpost <- pf_result$log_final_likelihood + lprior(prop_param) + ljacobian(prop_param);
  if(is.infinite(prop_lpost)) warning("initialization of the parameters get 0 likelihood!");
  prop_xT <- pf_result$xT;
  ## accept and reject
  if (log(runif(1)) < prop_lpost - curr_lpost){
    accept_chain[imcmc] <- 1;
    curr_param <- prop_param;
    curr_lpost <- prop_lpost;
    curr_xT <- prop_xT;
  }
  lpost_chain[imcmc] <- curr_lpost;
  ## save the parameters
  param_chain[imcmc, ] <- curr_param;
  xT_chain[imcmc, ] <- curr_xT;
  if(imcmc %% 100 ==0){
    cat("[",imcmc, "out of", num_mcmc, "iterations]\n");
    save.image("figures/section3/data_run_pmcmc_prediction_30obs_dgpinit.RData");
  }
}

cat("[ acceptance rate is ", mean(accept_chain), "]\n");
save.image("figures/section3/data_run_pmcmc_prediction_30obs_dgpinit.RData");
