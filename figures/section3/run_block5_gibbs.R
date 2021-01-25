rm(list = ls())
library(agents)
library(tictoc)
load("figures/section3/data_setup_hetero.RData");
source("figures/section3/setup_sis_pmcmc.R");


burn_in <- 5000 ;
num_mcmc <- 100 * 10**3 + burn_in;
dgp_parameters <- c(log(1/(N-1)), 0, dgp_bl1, dgp_bl2, dgp_bg1, dgp_bg2, dgp_config$rho);
step_size <- 0.08;
init_config <- dgp_config;

# expit <- function(z)  1 / (1 + exp(-z));
# get_lambda <- function(b1, b2) expit(b1 * features[1,] + b2 * features[2,]);

if(!is.null(init_config$policy)) init_config$policy <- NULL; ## gibbs does not need policy
### save and output these values
chain <- matrix(nrow = num_mcmc, ncol = 7); ## each row is a parameter
accept_chain <- rep(0, num_mcmc);
lpost_chain <- rep(NA, num_mcmc);
### instantiate parameters
set.seed(2020);
# curr_param <- c(rnorm(6,  sd = 3), runif(1)); ### RANDOM INIT
# curr_config <- update_model_config(init_config, curr_param); ### RANDOM INIT
# curr_param <- dgp_parameters; ### THIS INITIALIZES FROM DGP
# curr_config <- update_model_config(init_config, curr_param); ## this should be dgp_config
curr_param <- c(-10,4,-6,6,6,6,0.95); ## a bad point
curr_config <- update_model_config(init_config, curr_param);
init_param <- curr_param;
cat("initialized from", init_param, "\n");
curr_xx <- sis_xx_initialize(y, init_config);
curr_lpost <- sis_loglikelihood_complete(y, curr_xx, curr_config) + lprior(curr_param) + ljacobian(curr_param);
if(is.infinite(curr_lpost)) cat("initialization of the parameters get 0 likelihood!");
block_updates_state_space <- sis_get_state_space(5);

tic();
for(imcmc in 1 : num_mcmc){
  ## propose theta
  prop_param <- qkernel(curr_param, step_size);
  prop_config <- update_model_config(curr_config, prop_param);
  prop_lpost <- sis_loglikelihood_complete(y, curr_xx, prop_config) + lprior(prop_param) + ljacobian(prop_param);
  ## accept and reject theta
  if (log(runif(1)) < (prop_lpost - curr_lpost)){
    accept_chain[imcmc] <- 1
    curr_param <- prop_param;
    curr_lpost <- prop_lpost;
    curr_config <- prop_config;
  }
  ## update hidden states
  if(runif(1) < 0.5){
    ## swaps
    sis_xx_gibbs_swap_full_cpp(curr_xx, curr_config$alpha0, curr_config$lambda, curr_config$gamma);
  }else{
    ## block updates
    sis_xx_gibbs_blocked_full_cpp(curr_xx, y, curr_config$alpha0, curr_config$lambda, curr_config$gamma, curr_config$rho, 5, block_updates_state_space);
  }
  ## save the parameters
  chain[imcmc, ] <- curr_param;
  curr_lpost <- sis_loglikelihood_complete(y, curr_xx, curr_config) + lprior(curr_param) + ljacobian(curr_param);
  if(is.infinite(curr_lpost)) cat("blocked updates problem\n");
  lpost_chain[imcmc] <- curr_lpost;
  if(imcmc %% 1000 ==0){
    cat("[",imcmc, "out of", num_mcmc, "iterations]\n");
    save.image("figures/section3/data_run_block_gibbs_badinit.RData");
  } 
}
toc();
save.image("figures/section3/data_run_block_gibbs_badinit.RData");

# load("figures/section3/data_run_block_gibbs_test.RData")
# plot(lpost_chain, type = "l");
# cat("[acceptance probabilty = ", mean(accept_chain[1:imcmc]), "]\n");
# 
# plot.new()
# matplot(chain[10000: (imcmc - 1), ], col = 1: 7, type = "l", lty = 1)
# legend("topleft", paste(c(1:7)), col= c(1:7) , cex= 0.8,fill= c(1:7))
# abline(h = dgp_parameters, col =  1:7);
# 
# colMeans(chain[10000:(imcmc - 1), ])
# print(dgp_parameters)
