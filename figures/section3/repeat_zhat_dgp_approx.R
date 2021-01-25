rm(list = ls());
load("figures/section3/data_setup_hetero.RData");
library(agents)

## compute the backward information filter
dgp_config$policy <- sis_backward_information_filter_sumbin(y,dgp_config);

## setup particle filter
approx_config <- list(num_particles = 10,
                  ess_threshold = 0.5,
                  save_particles = FALSE,
                  clock = TRUE,
                  save_genealogy = FALSE,
                  verbose = FALSE,
                  exact = FALSE,
                  num_mcmc = ceiling(0.125 * dgp_config$N * log(dgp_config$N)));

# for different values of P, repeat the particle filters num_repeat times
num_particles_list <- 2**c(6,7,8,9,10,11);
num_repeats <- 100;
# num_particles_list <- 2**c(2,3); ## DEBUG
# num_repeats <- 10; ## DEBUG

df_length <- length(num_particles_list) * num_repeats * 2;
approx_llik_df_dgp <- data.frame(method = factor(character(df_length), levels = c("csmcApprox","apfApprox")),
                      nparticles = integer(df_length),
                      irep = integer(df_length),
                      logz = double(df_length),
                      time = double(df_length));

irow <- 1
for(p in num_particles_list){
  approx_config$num_particles <- p;
  for(irep in 1 : num_repeats){
    result_ <- sis_csmc(y, dgp_config, approx_config);
    approx_llik_df_dgp[irow,1] <- "csmcApprox";
    approx_llik_df_dgp[irow,2] <- p;
    approx_llik_df_dgp[irow,3] <- irep;
    approx_llik_df_dgp[irow,4] <- result_$log_final_likelihood;
    approx_llik_df_dgp[irow,5] <- result_$totaltime;
    irow <- irow + 1;
    cat("csmc repetition:", irep, "/", num_repeats, "particles:", p, "\n")
    
    result_ <- sis_apf(y, dgp_config, approx_config);
    approx_llik_df_dgp[irow,1] <- "apfApprox";
    approx_llik_df_dgp[irow,2] <- p;
    approx_llik_df_dgp[irow,3] <- irep;
    approx_llik_df_dgp[irow,4] <- result_$log_final_likelihood;
    approx_llik_df_dgp[irow,5] <- result_$totaltime;
    irow <- irow + 1;
    cat("apf repetition:", irep, "/", num_repeats, "particles:", p, "\n")
    
    # if(irow %% 100 ==0) cat("[",irow, "out of", df_length, "experiments done]\n");
  }
  save.image("figures/section3/data_repeat_zhat_dgp_approx.RData")
}

# head(approx_llik_df_dgp);
# tail(approx_llik_df_dgp);

