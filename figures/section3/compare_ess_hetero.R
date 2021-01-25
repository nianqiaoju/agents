rm(list = ls());
library(agents);
set.seed(2020); ##particle filters introduce randomness

load("figures/section3/data_setup_hetero.RData");

### some values for P, the number of particles
list_of_P <- 2 ** c(6:10);

### set up configuration for particle filters 
particle_config <- list(num_particles = 10,
                        save_particles = FALSE,
                        clock = FALSE,
                        save_genealogy = FALSE,
                        verbose = FALSE,
                        exact = TRUE);

### 3 methods: bpf, apf and csmc
### bif for csmc:
bif <- sis_backward_information_filter_sumbin(y, model_config = dgp_config);
model_config <- dgp_config;
model_config$policy <- bif;

### create a dataframe to store ess information
### NOTE::this is not good for memory but good enoughf for now. I NEED to optimize it later.
ess_df <- data.frame(method = factor(), P = integer(), step = integer(), ess = double());

for(np in list_of_P){
  ## update particle_config
  particle_config$num_particles <- np;
  ## bpf 
  bpf_ess <- sis_bpf(y,model_config,particle_config)$ess;
  ess_df <- rbind(ess_df, data.frame(method = factor("bpf"), P = np, step = 0 : (length(y) - 1), ess = bpf_ess));
  ## apf 
  apf_ess <- sis_apf(y,model_config,particle_config)$ess;
  ess_df <- rbind(ess_df, data.frame(method = factor("apf"), P = np, step = 0 : (length(y) - 1), ess = apf_ess));
  ## csmc
  csmc_ess <- sis_csmc_popobservation(y,model_config,particle_config)$ess;
  ess_df <- rbind(ess_df, data.frame(method = factor("csmc"), P = np, step = 0 : (length(y) - 1), ess = csmc_ess));
}


## quick ggplot
library(ggplot2)
unbiasedmcmc::setmytheme()
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
g <- ggplot(ess_df, aes(x = step, y = ess, colour = method)) + geom_line() + facet_grid(cols = vars(P));
g <- g + scale_colour_manual(values=cbPalette)
g
save.image("figures/section3/data_compare_ess_hetero.RData")
