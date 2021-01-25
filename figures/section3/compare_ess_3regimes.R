rm(list = ls());
library(agents);
library(tictoc)
set.seed(2020); ##particle filters introduce randomness

load("figures/section3/data_setup_hetero.RData");

### some values for P, the number of particles
num_particles <- 2 ** 9 
### set up configuration for particle filters 
particle_config <- list(num_particles = num_particles,
                        ess_threshold = 1,
                        save_particles = FALSE,
                        clock = FALSE,
                        save_genealogy = FALSE,
                        verbose = FALSE,
                        exact = TRUE);

### 3 regimes:
### original observations
### half the observations at time_steps
### double the observations at time_steps
mutated_y <- matrix(y, ncol =  3, nrow = length(y));
time_steps <- 1 + c(25,50,75);
mutated_y[time_steps, 2] <- floor(0.5 * y[time_steps]);
mutated_y[time_steps, 3] <- pmin(N, 2 * y[time_steps]);

colnames(mutated_y) <- c("original", "small", "large");
print(mutated_y[time_steps,]);

### create a dataframe to store ess information
### NOTE::this is not good for memory but good enoughffor now. I NEED to optimize it later.
ess_df <- data.frame(method = factor(), regime = factor(), step = integer(), ess = double());
model_config <- dgp_config;

for(regime in 1 : 3){
  cat("regime no.",regime,"\n");
  ### 4 methods: bpf, apf, csmc1 and csmc2

  ## bpf 
  bpf_ess <- sis_bpf(mutated_y[,regime],model_config,particle_config)$ess;
  ess_df <- rbind(ess_df, data.frame(method = factor("BPF"), regime = factor(regime), step = 0 : (length(y) - 1), ess = bpf_ess));
  
  ## apf 
  apf_ess <- sis_apf(mutated_y[,regime],model_config,particle_config)$ess;
  ess_df <- rbind(ess_df, data.frame(method = factor("APF"), regime = factor(regime), step = 0 : (length(y) - 1), ess = apf_ess));
  
  ## csmc with sumbin
  tic()
    bif <- sis_backward_information_filter_sumbin(mutated_y[,regime], model_config = dgp_config);
  timer <- toc()
  time_elapsed <- timer$toc - timer$tic
  cat("Runtime of BIF with SumBin:", time_elapsed, "\n")
  model_config$policy <- bif;
  csmc1_ess <- sis_csmc(mutated_y[,regime],model_config,particle_config)$ess;
  ess_df <- rbind(ess_df, data.frame(method = factor("cSMC1"), regime = factor(regime), step = 0 : (length(y) - 1), ess = csmc1_ess));
  
  ## csmc with transpoi
  tic()
    bif <- sis_backward_information_filter_tp(mutated_y[,regime], model_config = dgp_config);
  timer <- toc()
  time_elapsed <- timer$toc - timer$tic
  cat("Runtime of BIF with TransPoi:", time_elapsed, "\n")
  model_config$policy <- bif;
  csmc2_ess <- sis_csmc_tp(mutated_y[,regime],model_config,particle_config)$ess;
  ess_df <- rbind(ess_df, data.frame(method = factor("cSMC2"), regime = factor(regime), step = 0 : (length(y) - 1), ess = csmc2_ess));
  
}

save.image("figures/section3/data_compare_ess_3regimes.RData")
