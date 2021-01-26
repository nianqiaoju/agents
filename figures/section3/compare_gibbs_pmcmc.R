library(ggplot2);
library(dplyr);
rm(list = ls());
source("inst/set_plots.R"); ## set theme and colors

## compare pmmh, block5 and single-site


## load data ----
load("figures/section3/data_run_singlesite_gibbs_randominit.RData");
singlesite_rand_param_chain <- chain;
rm(chain);

load("figures/section3/data_run_block5_gibbs_randominit.RData");
block_rand_param_chain <- chain;
rm(chain);

load("figures/section3/data_run_pmcmc_randominitial.RData");
pmcmc_rand_param_chain <- chain;
rm(chain);

load("figures/section3/data_run_singlesite_gibbs_dgpinit.RData");
singlesite_dgp_param_chain <- chain;
rm(chain);

load("figures/section3/data_run_block5_gibbs_dgpinit.RData");
block_dpg_param_chain <- chain;
rm(chain);

load("figures/section3/data_run_pmcmc_final.RData");
pmcmc_dgp_param_chain <- chain;
rm(chain);

nburn <- burn_in;
nmcmc <- num_mcmc;

## compare the marginal distributions of beta[0]^1 and rho 
sample_index <- nburn : nmcmc;
sample_df <- data.frame(method = factor(character(), levels = c("pmmh", "block5","single-site")),
                        init = factor(character(), levels = c("DGP", "prior")),
                        index = integer(),
                        parameter_name = integer(),
                        parameter_value = double()
                        );


for(iparameter in c(1,7)){
  sample_df <- rbind(sample_df, data.frame(method = "pmmh", init = "prior",
                                           index = sample_index,
                                           parameter_name = iparameter,
                                           parameter_value = pmcmc_rand_param_chain[sample_index, iparameter]));
  sample_df <- rbind(sample_df, data.frame(method = "block5", init = "prior",
                                           index = sample_index,
                                           parameter_name = iparameter,
                                           parameter_value = block_rand_param_chain[sample_index, iparameter]));
  sample_df <- rbind(sample_df, data.frame(method = "single-site", init = "prior",
                                           index = sample_index,
                                           parameter_name = iparameter,
                                           parameter_value = singlesite_rand_param_chain[sample_index, iparameter]));
  sample_df <- rbind(sample_df, data.frame(method = "pmmh", init = "DGP",
                                           index = sample_index,
                                           parameter_name = iparameter,
                                           parameter_value = pmcmc_dgp_param_chain[sample_index, iparameter]));
  sample_df <- rbind(sample_df, data.frame(method = "block5", init = "DGP",
                                           index = sample_index,
                                           parameter_name = iparameter,
                                           parameter_value = block_dpg_param_chain[sample_index, iparameter]));
  sample_df <- rbind(sample_df, data.frame(method = "single-site", init = "DGP",
                                           index = sample_index,
                                           parameter_name = iparameter,
                                           parameter_value = singlesite_dgp_param_chain[sample_index, iparameter]));
}
head(sample_df);
summary(sample_df);

sample_df <- sample_df %>% mutate(parameter_name = factor(parameter_name, labels = c(expression(beta[0]^1), expression(rho))));

## create two histograms and paste them together
b0_hist <- ggplot(subset(sample_df, parameter_name == "beta[0]^1") , aes(x = parameter_value)) + geom_histogram(aes(y = stat(count) / sum(count)), bins = 12, fill = "grey60") 
b0_hist <- b0_hist + geom_vline(xintercept = dgp_parameters[1], size = 1) +   facet_grid( method ~ init , labeller = "label_parsed");
b0_hist <- b0_hist + ylab("density") + xlab(expression(beta[0]^1)) + ylim(0,0.06) +  theme(strip.text.y = element_blank());
# b0_hist;


rho_hist <- ggplot(subset(sample_df, parameter_name == "rho") , aes(x = parameter_value)) + geom_histogram(aes(y = stat(count) / sum(count)), bins = 15, fill = "grey60") 
rho_hist <- rho_hist + geom_vline(xintercept = dgp_parameters[7], size = 1) +   facet_grid( method ~ init , labeller = "label_parsed");
rho_hist <- rho_hist + xlab(expression(rho))	  + ylab(NULL) + ylim(0,0.06) + theme(axis.text.y = element_blank(), axis.title.x = element_text(margin=margin(32,0,0,0)));
# rho_hist;

library("cowplot");
p_hist1 <- plot_grid(b0_hist, rho_hist, ncol = 2, nrow = 1);


vline_df <- data.frame(parameter_name = levels(sample_df$parameter_name), value = dgp_parameters[c(1,7)]);
p_hist2 <- ggplot(sample_df , aes(x = parameter_value)) + geom_histogram(aes(y = stat(count) / sum(count)), bins = 12, fill = "grey60"); 
p_hist2 <- p_hist2 + facet_grid( method ~  parameter_name + init, labeller = "label_parsed", scales = "free_x");
p_hist2 <- p_hist2 + geom_vline(data = vline_df, aes(xintercept = value));
p_hist2 <- p_hist2 + ylab("density") + theme(axis.text.x = element_text(size = 14), axis.title.x = element_blank());
p_hist2

# ggsave(filename = "~/Dropbox/AgentBasedModels/paper/gibbs.compare.pmcmc.hist1.pdf", plot = p_hist1,
#        device = "pdf", width = 12, height = 12);
# ggsave(filename = "~/Dropbox/AgentBasedModels/paper/gibbs.compare.pmcmc.hist2.pdf", plot = p_hist2,
#        device = "pdf", width = 12, height = 12);
