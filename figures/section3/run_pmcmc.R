rm(list = ls())
library(agents)
library(tictoc)
load("figures/section3/data_setup_hetero.RData");
source("figures/section3/setup_sis_pmcmc.R");

burn_in <- 5000 
num_mcmc <- 100 * 10**3 + burn_in;
dgp_parameters <- c(log(1/(N-1)), 0, dgp_bl1, dgp_bl2, dgp_bg1, dgp_bg2, dgp_config$rho);

particle_config <- list(ess_threshold = 0.5,
                        num_particles = 128,
                        clock = FALSE,
                        save_genealogy = FALSE,
                        save_particles = FALSE,
                        verbose = FALSE,
                        exact = TRUE)
dgp_config$policy <- sis_backward_information_filter_sumbin(y, dgp_config);
step_size <- 0.2;
pf <- csmc;
init_config <- dgp_config;

# result <- pmcmc(num_mcmc, csmc, dgp_config, step_size);

### pf is a function that takes model_config as an input and returns logz 
### save and output these values
chain <- matrix(nrow = num_mcmc, ncol = 7); ## each row is a parameter
accept_chain <- rep(0, num_mcmc);
lpost_chain <- rep(0, num_mcmc);
### instantiate parameters
set.seed(2020);
curr_param <- c(rnorm(6,  sd = 3), runif(1)); ### RANDOM INIT
curr_config <- update_model_config(init_config, curr_param); ### RANDOM INIT
# curr_param <- dgp_parameters; ### THIS INITIALIZES FROM DGP
# curr_config <- dgp_config;
init_param <- curr_param;
print(init_param);
curr_lpost <- pf(curr_config) + lprior(curr_param) + ljacobian(curr_param);
if(is.infinite(curr_lpost)) stop("initialization of the parameters get 0 likelihood!");
# prop_config <- curr_config;
for(imcmc in 99101 : num_mcmc){
  ## propose
  prop_param <- qkernel(curr_param, step_size);
  prop_config <- update_model_config(curr_config, prop_param);
  prop_lpost <- pf(prop_config) + lprior(prop_param) + ljacobian(prop_param);
  ## accept and reject
  if (log(runif(1)) < prop_lpost - curr_lpost){
    accept_chain[imcmc] <- 1;
    curr_param <- prop_param;
    curr_lpost <- prop_lpost;
  }
  lpost_chain[imcmc] <- curr_lpost;
  ## save the parameters
  chain[imcmc, ] <- curr_param;
  if(imcmc %% 100 ==0){
    cat("[",imcmc, "out of", num_mcmc, "iterations]\n");
    save.image("figures/section3/data_run_pmcmc_randominitial.RData");
  }
}
save.image("figures/section3/data_run_pmcmc_randominitial.RData");

# plot(lpost_chain, type = "l");

# ### 
# plot.new()
# matplot(chain, col = 1: 7, type = "l", lty = 1)
# legend("topleft", paste(c(1:7)), col= c(1:7) , cex= 0.8,fill= c(1:7))
# abline(h = dgp_parameters, col =  1:7);
# mean(accept_chain)
# 
# colMeans(chain[5000:(imcmc - 1), ])
# print(dgp_parameters)
