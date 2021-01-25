## compute the variance of log-marginal likelihood estimator and average cost
## for different values of P and SMC methods(apf, bpf, csmc1, csmc2)
## at dgp and a non-dgp value
## the table "variance_summary" becomes Table 2 in Section 3.3

library(ggplot2)
library(dplyr)

# exact SMC methods at DGP
rm(list = ls())
load("figures/section3/data_repeat_zhat_dgp_exact.RData") 
cost_p <- ggplot(data = llik_df, aes(x = time, group = method)) + geom_histogram() + facet_grid( . ~ method)
cost_p

logzhat_p <- ggplot(data = llik_df, aes(x = logz)) + geom_histogram() + facet_grid(nparticles ~ method)
logzhat_p

variance_summary <- llik_df %>% group_by(method, nparticles) %>% summarise(variance = var(logz), cost = mean(time))
variance_summary %>% mutate(inefficiency = variance * cost) -> variance_summary
variance_summary %>% mutate(asymptoticvar = nparticles * variance) -> variance_summary
variance_summary

write.table(variance_summary, file = "figures/section3/results_zhat_dgp_exact.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)


# exact SMC methods at non-DGP
rm(list = ls())
load("figures/section3/data_repeat_zhat_nondgp_exact.RData") 
cost_p <- ggplot(data = llik_df, aes(x = time, group = method)) + geom_histogram() + facet_grid( . ~ method)
cost_p

logzhat_p <- ggplot(data = llik_df, aes(x = logz)) + geom_histogram() + facet_grid(nparticles ~ method)
logzhat_p

variance_summary <- llik_df %>% group_by(method, nparticles) %>% summarise(variance = var(logz), cost = mean(time))
variance_summary %>% mutate(inefficiency = variance * cost) -> variance_summary
variance_summary %>% mutate(asymptoticvar = nparticles * variance) -> variance_summary
variance_summary

write.table(variance_summary, file = "figures/section3/results_zhat_nondgp_exact.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE)



