rm(list = ls());
load("figures/section3/data_setup_hetero.RData");

### create a grid for rho
rho_grid <- seq(0.1, 0.9, length.out = 9);
num_repeats <- 100;
num_particles <- 2^9 

### 
csmc_config  <- list(num_particles = num_particles,
                     ess_threshold = 0.5,
                     save_particles = FALSE,
                     clock = FALSE,
                     save_genealogy = FALSE,
                     verbose = FALSE,
                     exact = TRUE);

model_config <- dgp_config;

df_length <- length(rho_grid) * num_repeats * 4;
grid_rho_df <- data.frame(rho = double(df_length),
                          method = factor(character(df_length), levels = c("csmc", "csmc-tp", "apf", "bpf")),
                          irep = integer(df_length),
                          llik = double(df_length));
irow <- 1;
for(rho in rho_grid){
  model_config$rho <- rho;
  model_config$policy <- sis_backward_information_filter_sumbin(y,  model_config);
  for(irep in 1 : num_repeats){
    zhat <- sis_csmc(y, model_config, csmc_config)$log_final_likelihood;
    grid_rho_df[irow,1] <- rho;
    grid_rho_df[irow,2] <- "csmc";
    grid_rho_df[irow,3] <- irep;
    grid_rho_df[irow,4] <- zhat;
    irow <- irow + 1;
    cat("csmc repetition:", irep, "/", num_repeats, "rho:", rho, "\n")
    cat("log-likelihood estimate:", zhat, "\n")
    
    zhat <- sis_csmc_tp(y, model_config, csmc_config)$log_final_likelihood;
    grid_rho_df[irow,1] <- rho;
    grid_rho_df[irow,2] <- "csmc-tp";
    grid_rho_df[irow,3] <- irep;
    grid_rho_df[irow,4] <- zhat;
    irow <- irow + 1;
    cat("csmc-tp repetition:", irep, "/", num_repeats, "rho:", rho, "\n")
    cat("log-likelihood estimate:", zhat, "\n")
    
    zhat <- sis_bpf(y, model_config, csmc_config)$log_final_likelihood;
    grid_rho_df[irow,2] <- "bpf";
    grid_rho_df[irow,1] <- rho;
    grid_rho_df[irow,3] <- irep;
    grid_rho_df[irow,4] <- zhat
    irow <- irow + 1;
    cat("bpf repetition:", irep, "/", num_repeats, "rho:", rho, "\n")
    cat("log-likelihood estimate:", zhat, "\n")
    
    zhat <- sis_apf(y, model_config, csmc_config)$log_final_likelihood;
    grid_rho_df[irow,2] <- "apf";
    grid_rho_df[irow,1] <- rho;
    grid_rho_df[irow,3] <- irep;
    grid_rho_df[irow,4] <- zhat
    irow <- irow + 1;
    cat("apf repetition:", irep, "/", num_repeats, "rho:", rho, "\n")
    cat("log-likelihood estimate:", zhat, "\n")
  }
  save.image("figures/section3/data_generate_grid_rho_fix_beta.RData");
}
