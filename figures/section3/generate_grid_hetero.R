rm(list = ls());
load("figures/section3/data_setup_hetero.RData");


### create a grid for beta_lambda_1 & beta_lambda_2 and check for identifiablity
bl1_grid <- bl2_grid <- seq(-6,6,length.out = 30);
time_grid <- c(10,30,90);

num_particles <- 2**6; ## 
### use csmc until the final choice is made
csmc_config  <- list(num_particles = num_particles,
                     ess_threshold = 0.5,
                    save_particles = FALSE,
                    clock = FALSE,
                    save_genealogy = FALSE,
                    verbose = FALSE,
                    exact = TRUE);

## create an empty dataframe of the right size
df_length <- length(time_grid) * length(bl1_grid) * length(bl2_grid);

llik_df <- data.frame(t = integer(df_length), 
                      bl1 = double(df_length), 
                      bl2 = double(df_length),
                      llik = double(df_length));

model_config <- dgp_config;

irow <- 1;
for(bl1 in bl1_grid){
  for(bl2 in bl2_grid){
    ## modify configuration
    model_config$lambda <- get_lambda(b1 = bl1, b2 = bl2)
    for(time_step in time_grid){
      ## compute bif 
      bif <- sis_backward_information_filter_sumbin(y[1 : (1 + time_step)], model_config);
      model_config$policy <- bif;
      ## compute logz
      logz <- sis_csmc(y = y[1 : (1 + time_step)], model_config  = model_config, particle_config = csmc_config)$log_final_likelihood;
      ## store information in llik_df
      llik_df[irow, 1] <- time_step;
      llik_df[irow, 2] <- bl1;
      llik_df[irow, 3] <- bl2;
      llik_df[irow, 4] <- logz;
      irow <- irow + 1;
      if(irow %% 100 ==0) cat("[",irow, "out of", df_length, "experiments done]\n");
    }
  }
}

save.image("figures/section3/data_generate_grid_hetero.RData")

