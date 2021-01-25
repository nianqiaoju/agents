# rm(list = ls());
# load("figures/section3/data_setup_hetero.RData");
# library(agents)
# 
# new_config <- dgp_config
# new_lambda <- get_lambda(-3,0);
# new_config$lambda <- new_lambda
# 
# ## compute the backward information filter
# new_config$policy <- sis_backward_information_filter_sumbin(y,new_config);
# 
# ## setup particle filter
# pf_config <- list(num_particles = 10,
#                   ess_threshold = 0.5,
#                   save_particles = FALSE,
#                   clock = TRUE,
#                   save_genealogy = FALSE,
#                   verbose = FALSE, 
#                   exact = TRUE)
# 
# # for different values of P, repeat the particle filters num_repeat times
# num_particles_list <- 2**c(6,7,8,9,10,11);
# num_repeats <- 100;
# # num_particles_list <- 2**c(2,3); ## DEBUG
# # num_repeats <- 1; ## DEBUG
# 
# df_length <- length(num_particles_list) * num_repeats * 3;
# llik_df <- data.frame(method = factor(character(df_length), levels = c("csmc","apf","bpf")),
#                       nparticles = integer(df_length),
#                       irep = integer(df_length),
#                       logz = double(df_length),
#                       time = double(df_length));
# 
# irow <- 1
# for(p in num_particles_list){
#   pf_config$num_particles <- p;
#   for(irep in 1 : num_repeats){
#     result_ <- sis_csmc(y, new_config, pf_config);
#     llik_df[irow,1] <- "csmc";
#     llik_df[irow,2] <- p;
#     llik_df[irow,3] <- irep;
#     llik_df[irow,4] <- result_$log_final_likelihood;
#     llik_df[irow,5] <- result_$totaltime;
#     irow <- irow + 1;
#     cat("csmc repetition:", irep, "/", num_repeats, "particles:", p, "\n")
#     
#     result_ <- sis_bpf(y, new_config, pf_config);
#     llik_df[irow,1] <- "bpf";
#     llik_df[irow,2] <- p;
#     llik_df[irow,3] <- irep;
#     llik_df[irow,4] <- result_$log_final_likelihood;
#     llik_df[irow,5] <- result_$totaltime;
#     irow <- irow + 1;
#     cat("bpf repetition:", irep, "/", num_repeats, "particles:", p, "\n")
#     
#     result_ <- sis_apf(y, new_config, pf_config);
#     llik_df[irow,1] <- "apf";
#     llik_df[irow,2] <- p;
#     llik_df[irow,3] <- irep;
#     llik_df[irow,4] <- result_$log_final_likelihood;
#     llik_df[irow,5] <- result_$totaltime;
#     irow <- irow + 1;
#     cat("apf repetition:", irep, "/", num_repeats, "particles:", p, "\n")
#     
#     # if(irow %% 100 ==0) cat("[",irow, "out of", df_length, "experiments done]\n");
#   }
#   save.image("figures/section3/data_repeat_zhat_nondgp_exact.RData")
# }
# 
# # head(llik_df);
# # tail(llik_df);

# extra repetitions for bpf with increased no. of particles 
rm(list = ls());
load("figures/section3/data_repeat_zhat_nondgp_exact.RData")
bpf_llike <- rep(0, num_repeats)
bpf_time <- rep(0, num_repeats)
pf_config$num_particles <- 2^18
for (irep in 1 : num_repeats){
  result_ <- sis_bpf(y, new_config, pf_config)
  bpf_llike[irep] <- result_$log_final_likelihood
  bpf_time[irep] <- result_$totaltime
  if (irep %% 10 == 0){
    save.image("figures/section3/data_repeat_zhat_nondgp_exact.RData")
  }
  cat("bpf repetition:", irep, "/", num_repeats, "\n")
}