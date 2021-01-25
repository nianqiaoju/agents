## this script runs a long MH-exact chain
## for the set up of N = 1000
rm(list = ls())
library(agents)
load(file = "figures/section2/data_tune_mh.RData")
set.seed(2020)

## mcmc setup 
num_burn <- 5000
num_mcmc_long <- 100000 + num_burn
step_size <- 0.2

## features of each agent
n <- 1000;
features <- config_list[[2]]$features
model_config <- list(features = features)
## observation
y <- y[2];

## define functions to use for likelihood etc. 
llik_exact <- function(parameters){
  static_loglikelihood(y = y, parameters = parameters, model_config = model_config, method = "exact")
}

cat("[ running long exact mh chain for population size", n, "]\n")
## exact mh
result_ <- mh(num_mcmc_long, llik_exact)
beta_star <- mean(result_$beta_chain[num_burn : num_mcmc_long]);
rho_star <- mean(result_$rho_chain[num_burn : num_mcmc_long]);

save.image(file = "figures/section2/data_posterior_long_n1000.RData");
