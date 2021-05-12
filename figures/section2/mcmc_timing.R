rm(list = ls())
library(agents)
load(file = "figures/section2/data_tune_mh.RData")
set.seed(2020)

num_repeats <- 1
## mcmc setup 
num_burn <- 500
num_mcmc_short <- 500 + num_burn
# num_mcmc_long <- 50000 + num_burn
step_size <- 0.2

#### repeated short runs for each data config
post_df <- data.frame()
for (irep in 1 : num_repeats){
  cat("[running repeat no.", irep,"]\n")
  for (i in 1 : size_N){
    cat("[population size", sub_N[i],"]\n")
    
    ## define functions to use for likelihood etc. 
    llik_exact <- function(parameters){
      static_loglikelihood_marginal(y = y[i], parameters = parameters, model_config = config_list[[i]], method = "exact")
    }
    llik_tp <- function(parameters){
      static_loglikelihood_marginal(y = y[i], parameters = parameters, model_config = config_list[[i]], method = "tp")
    }
    llik_is <- function(parameters){
      static_loglikelihood_marginal(y = y[i], parameters = parameters, model_config = config_list[[i]], method = "mc", particle_config = particle_config)
    }
    llik_complete <- function(parameters, xx){
      alpha_ <- get_rates_from_features(beta = parameters$beta, features = config_list[[i]]$features)
      dbinom(x = y[i], size = sum(xx), prob = parameters$rho, log = T) + lw.logsum(dbinom(x = xx, size = 1, prob = alpha_, log = T))
    }
    xx_gibbs_scan <- function(parameters,xx){
      rho_ <- parameters$rho
      alpha_ <- get_rates_from_features(beta = parameters$beta, features = config_list[[i]]$features)
      return(as.numeric(static_xx_gibbs_cpp(xx_previous = xx, alpha = alpha_, rho = parameters$rho, y = y[i])))
    }
    xx_gibbs_block <- function(parameters, xx){
      rho_ <- parameters$rho
      alpha_ <- get_rates_from_features(beta = parameters$beta, features = config_list[[i]]$features)
      return(as.numeric(static_xx_blocked_gibbs(xx = xx, alpha = alpha_, rho = parameters$rho, y = y[i], block_size = 20)))
    }
    xx_swap <- function(parameters, xx){
      alpha_<- get_rates_from_features(beta = parameters$beta, features = config_list[[i]]$features)
      return(as.numeric(static_xx_swap(xx, alpha_)))
    }
    ## exact mh
    result1_ <- mh(num_mcmc_short, llik_exact)
    post_df <- rbind(post_df, data.frame(N = sub_N[i], 
                                         method = factor("mh-exact"),
                                         Ebeta = mean(result1_$beta_chain[num_burn : num_mcmc_short]),
                                         Erho = mean(result1_$rho_chain[num_burn : num_mcmc_short]),
                                         correlation = cor(result1_$beta_chain[num_burn : num_mcmc_short], result1_$rho_chain[num_burn : num_mcmc_short]),
                                         cost = result1_$elapse))
    
    ## mh with biased likelihood
    result2_ <- mh(num_mcmc_short, llik_tp)
    post_df <- rbind(post_df, data.frame(N = sub_N[i], 
                                         method = factor("mh-tp"),
                                         Ebeta = mean(result2_$beta_chain[num_burn : num_mcmc_short]),
                                         Erho = mean(result2_$rho_chain[num_burn : num_mcmc_short]),
                                         correlation = cor(result2_$beta_chain[num_burn : num_mcmc_short], result2_$rho_chain[num_burn : num_mcmc_short]),
                                         cost = result2_$elapse))
    
    ## pmmh using importance sampling
    result3_ <- mh(num_mcmc_short, llik_is)
    post_df <- rbind(post_df, data.frame(N = sub_N[i], 
                                         method = factor("pmmh"),
                                         Ebeta = mean(result3_$beta_chain[num_burn : num_mcmc_short]),
                                         Erho = mean(result3_$rho_chain[num_burn : num_mcmc_short]),
                                         correlation = cor(result3_$beta_chain[num_burn : num_mcmc_short], result3_$rho_chain[num_burn : num_mcmc_short]),
                                         cost = result3_$elapse))
    
    ## single-site gibbs
    result4_ <- gibbs(num_mcmc_short, llik_complete, xx_gibbs_scan, xx_swap, sub_N[i])
    post_df <- rbind(post_df, data.frame(N = sub_N[i], 
                                         method = factor("gibbs-singlesite"),
                                         Ebeta = mean(result4_$beta_chain[num_burn : num_mcmc_short]),
                                         Erho = mean(result4_$rho_chain[num_burn : num_mcmc_short]),
                                         correlation = cor(result4_$beta_chain[num_burn : num_mcmc_short], result4_$rho_chain[num_burn : num_mcmc_short]),
                                         cost = result4_$elapse))
    
    ## blocked gibbs
    result5_ <- gibbs(num_mcmc_short, llik_complete, xx_gibbs_block, xx_swap, sub_N[i])
    post_df <- rbind(post_df, data.frame(N = sub_N[i], 
                                         method = factor("gibbs-blocked"),
                                         Ebeta = mean(result5_$beta_chain[num_burn : num_mcmc_short]),
                                         Erho = mean(result5_$rho_chain[num_burn : num_mcmc_short]),
                                         correlation = cor(result5_$beta_chain[num_burn : num_mcmc_short], result5_$rho_chain[num_burn : num_mcmc_short]),
                                         cost = result5_$elapse))
  }
}

head(post_df)
tail(post_df)
View(post_df)
post_df$cost <- post_df$cost / 10
save.image(file = "figures/section2/mcmc_timing.RData")
