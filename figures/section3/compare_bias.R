rm(list = ls())
library(ggplot2)
library(dplyr)

load("figures/section3/data_repeat_zhat_dgp_exact.RData") ## dgp exact
dgp_exact_llik <- llik_df;
load("figures/section3/data_repeat_zhat_dgp_approx.RData") # non-DGP
dgp_approx_llik <- approx_llik_df_dgp;

### unbiased estimator of z:
logzs <- dgp_exact_llik$logz;
maxlogzs <- max(logzs);
logzs <- logzs - maxlogzs;
zhat <- mean(exp(logzs)) * exp(maxlogzs);
cat("[an unbiased estimate of z is ", zhat, "]\n");

### helper function: 
## input: a vector of estimations of logz;
## output: bias / zhat
get_bias <- function(logz){
  maxlz <- max(logz);
  abs(mean(exp(logz - maxlz)) * exp(maxlz) - zhat) / zhat;
}


dgp_approximation_bias <- dgp_approx_llik %>% group_by(method, nparticles) %>% summarise(bias =  get_bias(logz));
dgp_approximation_bias

# write.table(variance_summary, file = "figures/section3/output_variance_logzhat.txt", sep = "\t",
            # row.names = FALSE, col.names = TRUE)



