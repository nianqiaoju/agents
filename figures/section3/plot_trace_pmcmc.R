rm(list = ls());
# load("figures/section3/data_run_pmcmc_test.RData");
# load("figures/section3/data_run_pmcmc_dgp_init.RData");
load("figures/section3/data_run_pmcmc_final.RData");
# trace plots
plot.new()
matplot(chain, col = 1: 7, type = "l", lty = 1)
legend("topleft", paste(c(1:7)), col= c(1:7) , cex= 0.8,fill= c(1:7))
abline(h = dgp_parameters, col =  1:7);

# average acceptance probability (seems acceptable to me, I think no need to tune proposal SD)
mean(accept_chain)

# log-posterior trace (chain is moving but quite sticky, might want to increase P?)
plot(lpost_chain[1:1000], type = "l") 
plot(lpost_chain[10000:11000], type = "l") 
plot(lpost_chain[20000:21000], type = "l") 
plot(lpost_chain[30000:31000], type = "l") 
plot(lpost_chain[40000:41000], type = "l") 

# posterior means with burn-in = 5000
colMeans(chain[(burn_in+1):(num_mcmc), ])

# DGP
print(dgp_parameters)

# posterior marginal distributions (the parameters that the means do not recover are quite uncertain)
hist(chain[,1], breaks = 50, probability = TRUE) # unimodal, left-skewed, mean close to DGP
hist(chain[,2], breaks = 50, probability = TRUE) # unimodal, quite symmetric, mean close to DGP
hist(chain[,3], breaks = 50, probability = TRUE) # unimodal, right-skewed, mean positive while DGP is negative
hist(chain[,4], breaks = 50, probability = TRUE) # unimodal, slightly right-skewed, mean less positive than DGP
hist(chain[,5], breaks = 50, probability = TRUE) # unimodal, slightly right-skewed, mean positive while DGP is negative
hist(chain[,6], breaks = 50, probability = TRUE) # multi-modal, mean suprisngly close to DGP  
hist(chain[,7], breaks = 50, probability = TRUE) # very right-skewed, mean close to DGP  

# auto-correlation plots with burn-in = 5000 (mixing is quite bad)
acf(chain[(burn_in+1):num_mcmc, 1], lag.max = 100)
acf(chain[(burn_in+1):num_mcmc, 2], lag.max = 100)
acf(chain[(burn_in+1):num_mcmc, 3], lag.max = 100)
acf(chain[(burn_in+1):num_mcmc, 4], lag.max = 100)
acf(chain[(burn_in+1):num_mcmc, 5], lag.max = 100)
acf(chain[(burn_in+1):num_mcmc, 6], lag.max = 100)
acf(chain[(burn_in+1):num_mcmc, 7], lag.max = 100)

