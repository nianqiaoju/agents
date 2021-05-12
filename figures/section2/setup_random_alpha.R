rm(list =ls())
library(agents)
set.seed(2020)

## parameters to simulate y 
N <- 1500 # population size
sub_N <- c(500,1000,1500)
size_N <- length(sub_N)
num_features <- 1
features <- matrix(rnorm(N, mean = 4), nrow = num_features) # random features
dgp <- list(beta = 0.3, rho = 0.8)

## compute alpha
alpha <- get_rates_from_features(features = features, beta =  dgp$beta)

## make model_config for each sub-population
config_list <- list()
for (i in 1 : length(sub_N)){
  config_list[[i]] <- list(N = sub_N[i], features = matrix(features[1, 1 : sub_N[i]], nrow = 1))
}

## simulate observations in sub populations
z <- runif(N)
y <- sapply(sub_N, function(n) sum(z[1:n] < alpha[1:n] * dgp$rho))
print(y)

rm(i)
save.image(file = 'figures/section2/data_static_observations.RData')
