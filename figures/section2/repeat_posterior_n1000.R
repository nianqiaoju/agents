## this script measures the monte carlo error of mh-exact, mh-tp, pmmh and gibbs-full
## it repeats short runs of each algorithm for 25 times
## and record the posterior mean of each parameter from each run

rm(list = ls())
library(agents)
source("figures/section2/tune_mh.R")
set.seed(2020)

num_repeats <-  25
## mcmc setup 
num_burn <- 5000
num_mcmc_short <- 20000 + num_burn
num_mcmc_long <- 50000 + num_burn
step_size <- 0.2

features <- config_list[[2]]$features
n <- 1000;
y <- y[2];
model_config <- list(features = features)

#### repeated short runs for each data config
post_df <- data.frame()
for (irep in 1 : num_repeats){
  cat("[running repeat no.", irep,"]\n")

  ## define functions to use for likelihood etc. 
  llik_exact <- function(parameters){
    static_loglikelihood(y = y, parameters = parameters, model_config = model_config, method = "exact")
  }
  llik_tp <- function(parameters){
    static_loglikelihood(y = y, parameters = parameters, model_config = model_config,  method = "tp")
  }
  llik_is <- function(parameters){
    static_loglikelihood(y = y, parameters = parameters, model_config = model_config, method = "mc", particle_config = particle_config)
  }
  llik_complete <- function(parameters, xx){
    alpha_ <- static_get_alpha(beta = parameters$beta, features = features)
    dbinom(x = y, size = sum(xx), prob = parameters$rho, log = T) + lw.logsum(dbinom(x = xx, size = 1, prob = alpha_, log = T))
  }
  xx_gibbs_scan <- function(parameters,xx){
    rho_ <- parameters$rho
    alpha_ <- static_get_alpha(beta = parameters$beta, features = features)
    return(as.numeric(static_xx_gibbs_cpp(xx_previous = xx, alpha = alpha_, rho = parameters$rho, y = y)))
  }
  xx_swap <- function(parameters, xx){
    alpha_<- static_get_alpha(beta = parameters$beta, features = features)
    return(as.numeric(static_xx_swap(xx, alpha_)))
  }
  ## exact mh
  result1_ <- mh(num_mcmc_short, llik_exact)
  post_df <- rbind(post_df, data.frame(N = n, 
                                       method = factor("mh-exact"),
                                       Ebeta = mean(result1_$beta_chain[num_burn : num_mcmc_short]),
                                       Erho = mean(result1_$rho_chain[num_burn : num_mcmc_short]),
                                       cost = result1_$elapse))
  
  ## mh with biased likelihood
  result2_ <- mh(num_mcmc_short, llik_tp)
  post_df <- rbind(post_df, data.frame(N = n, 
                                       method = factor("mh-tp"),
                                       Ebeta = mean(result2_$beta_chain[num_burn : num_mcmc_short]),
                                       Erho = mean(result2_$rho_chain[num_burn : num_mcmc_short]),
                                       cost = result2_$elapse))
  
  ## pmmh using importance sampling
  result3_ <- mh(num_mcmc_short, llik_is)
  post_df <- rbind(post_df, data.frame(N = n, 
                                       method = factor("pmmh"),
                                       Ebeta = mean(result3_$beta_chain[num_burn : num_mcmc_short]),
                                       Erho = mean(result3_$rho_chain[num_burn : num_mcmc_short]),
                                       cost = result3_$elapse))
  
  ## single-site gibbs
  result4_ <- gibbs(num_mcmc_short, llik_complete, xx_gibbs_scan, xx_swap, n)
  post_df <- rbind(post_df, data.frame(N = n, 
                                       method = factor("gibbs-singlesite"),
                                       Ebeta = mean(result4_$beta_chain[num_burn : num_mcmc_short]),
                                       Erho = mean(result4_$rho_chain[num_burn : num_mcmc_short]),
                                       cost = result4_$elapse))
  save.image(file = "figures/section2/data_repeat_posterior_n1000.RData")
}

head(post_df)
tail(post_df)
View(post_df)

save.image(file = "figures/section2/data_repeat_posterior_n1000.RData")

