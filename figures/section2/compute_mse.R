## this script computes the Monte Carlo error of estimating 
## the posterior mean of beta and rho 
## the table bv_df is Table 1 in the paper.
rm(list = ls());
# load("figures/section2/data_sample_posterior_repeats.RData")
# post_df <- subset(post_df, N == 1000);
load("figures/section2/data_repeat_posterior_n1000.RData");
load("figures/section2/data_posterior_long.RData");
library(dplyr);

### separate bias and variance
bv_df <- post_df %>% group_by(method) %>% summarise(beta_mean = mean(Ebeta),
                                                   beta_variance = var(Ebeta),
                                                   rho_mean = mean(Erho),
                                                   rho_variance = var(Erho),
                                                   cost = mean(cost));
## these are posterior mean from a long run of the MH-exact algorithm
cat("posterior mean of beta is", beta_star);
cat("posterior mean of rho is", rho_star);

bv_df <- bv_df %>% mutate(beta_bias2 = (beta_mean - beta_star)**2, rho_bias2 = (rho_mean - rho_star)**2);
bv_df <- bv_df %>% mutate(beta_mse = beta_variance + beta_bias2, rho_mse = rho_variance + rho_bias2);
bv_df$beta_mean <- NULL;
bv_df$rho_mean <- NULL;
bv_df

### auto correlation plot
acf(x = result1_$beta_chain, plot = F)
acf(x = result2_$beta_chain, plot = F)
acf(x = result3_$beta_chain, plot = F)
acf(x = result4_$beta_chain, plot = F)
# acf(x = result5_$beta_chain, plot = F)
