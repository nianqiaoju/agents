rm(list = ls())
set.seed(2020)
library(agents)
## set up the features and network
N <- 100;
steps <- 90;
network <- network_fully_connected(N);
## features for heterogenous agents
features <- matrix(nrow = 2, ncol = N);
features[1,] <- 1;
# features[2,] <- rgamma(N, shape = 3);
# features[2,] <- (features[2,] - mean(features[2,]))/sd(features[2,]);
features[2,] <- rnorm(N);

## make lambda and gamma depend on features
dgp_bl1 <- - 1;
dgp_bl2 <- 2; ## this is signal strength
get_lambda <- function(b1,b2){
  expit(b1 * features[1,] + b2 * features[2,]);
}
dgp_lambda <- get_lambda(dgp_bl1,dgp_bl2);
dgp_bg1 <- - 1;
dgp_bg2 <- - 1;
dgp_gamma <- get_lambda(dgp_bg1, dgp_bg2);

## make sure there is enough heterogenity in lambda and gamma
###
cat("[data generating lambda ranges from", min(dgp_lambda),"to", max(dgp_lambda),"]\n");
cat("[with mean", mean(dgp_lambda),"sd", sd(dgp_lambda),"]\n");
cat("[data generating gamma ranges from", min(dgp_gamma),"to", max(dgp_gamma),"]\n");
cat("[with mean", mean(dgp_gamma),"sd", sd(dgp_gamma),"]\n");
hist(dgp_lambda);
hist(dgp_gamma);

dgp_config <- list(N = N, 
                  adjacency_matrix_b = network,
                  network_type = 'full',
                   alpha0 = rep(1/ N,N),
                   gamma = dgp_gamma,
                   lambda = dgp_lambda,
                   rho = 0.8);

## simulate observations
complete_data <- sis_simulate(days = steps, model_config = dgp_config, dt = 1);
y <- complete_data$y;
## reject the process if observations are zero
while(sum(y <= 1) > 10){
  complete_data <- sis_simulate(days = steps, model_config = dgp_config, dt = 1)
  y <- y <- complete_data$y;
}

plot(y, type = "l");

## save
save.image("figures/section3/data_setup_hetero.RData");
