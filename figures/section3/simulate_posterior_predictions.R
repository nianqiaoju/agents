rm(list = ls());
library(ggplot2);
library(agents);
source("inst/set_plots.R")
load("figures/section3/data_run_pmcmc_prediction_30obs_dgpinit.RData");
set.seed(2020);

## predict number of observed patients on each day ----
sis_predict <- function(xt, curr_param, days = (length(y) - 1 - obs_length)){
  ypred <- rep(NA, days);
  lambda <- get_lambda(curr_param[3], curr_param[4]);
  gamma <- get_lambda(curr_param[5], curr_param[6]);
  alpha_t <- sis_get_alpha_full_cpp(xt, lambda, gamma);
  for(d in 1 : days){
    agent_states <- (runif(N) < alpha_t);
    ypred[d] <- rbinom(1, sum(agent_states), curr_param[7]);
    alpha_t <- sis_get_alpha_full_cpp(agent_states, lambda, gamma);
  }
  return(ypred)
}

nmcmc <- 105000
nburn <- 55000

ypred_samples <- matrix(NA, (length(y) - obs_length - 1), nmcmc - nburn);
for(isample in (nburn + 1) : nmcmc){
  ## simulate a trajectory
  ypred_samples[, isample - nburn] <- sis_predict(xt = xT_chain[isample, ], curr_param = param_chain[isample,]);
}


## distribution of predicted y[31] given y[0:30];
hist(ypred_samples[1,]);

## how many of the predictions leads to no breakout?
no_outbreaks <- function(ypred){
  sum(ypred == 0) >= 2;
}

## ~43 of the predicted trajectories are not outbreaks
outbreaks_id <- which(!apply(ypred_samples, 2 , no_outbreaks));
print(1 - length(outbreaks_id) / (nmcmc - nburn));


## a dataframe for the observed ys
data_df <- data.frame(x = 0 : (length(y) - 1 ), y = y, observed = c(rep(TRUE, obs_length + 1), rep(FALSE,  length(y) - obs_length - 1)));
data_df$observed <- factor(data_df$observed);

## a dataframe for the trajectories
traj_df <- data.frame(x = integer(length(y) - obs_length - 1), y = integer(length(y) - obs_length - 1));
traj_df$x <-  (obs_length  + 1 ) : (length(y) - 1 )
head(traj_df);
## thin the trajectories for visualization
sample_id <- as.integer(seq(from = 1, to = nmcmc - nburn, length.out = 100));


## get median and quantiles y for each day
get_quantiles <- function(xx){
  quantile(xx, probs = c(0.025,0.5, 1 - 0.025));
}

y_pred_quantiles <- apply(ypred_samples, 1, get_quantiles);
# ypred_nozero_quantiles <- apply(ypred_samples[, outbreaks_id], 1, get_nozero_quantile);

## a dataframe for the quantiles
quantile_df <- data.frame(x = (obs_length  + 1 ) : (length(y) - 1 ), 
                          ypredmin = y_pred_quantiles[1, ], 
                          ypredmed = y_pred_quantiles[2, ],
                          ypredmax = y_pred_quantiles[3, ]);

## prepare plot for paper-----
pred_plot <- ggplot(data_df, aes(x = x, y = y));
pred_plot <- pred_plot + geom_point(size = 3, shape = 21, aes(fill = observed));
## plot the predicted trajctories
for(isample in sample_id){
  traj_df$y <- ypred_samples[,  isample];
  pred_plot <- pred_plot + geom_line(data = traj_df, aes(x = x, y = y), color = "grey80");
}
## add the quantile lines
pred_plot <- pred_plot + geom_line(data = quantile_df, aes(x = x, y = ypredmed), linetype = 2, size = 2, color = "blue") + 
  geom_line(data = quantile_df, aes(x = x, y = ypredmin), linetype = 2, size = 2, color = cbPalette[3]) +
  geom_line(data = quantile_df, aes(x = x, y = ypredmax), linetype = 2, size = 2, color = cbPalette[3]);
pred_plot <- pred_plot + geom_line(data = data_df[(obs_length + 2) : length(y), ]);
## further formatting
pred_plot <- pred_plot + xlab("time") + ylim(0,80) + ylab("observation");
pred_plot <- pred_plot + geom_point(size = 3, shape = 21, aes(fill = observed));
pred_plot <- pred_plot + scale_fill_grey(start = 0.8, end = 0) + theme(legend.position = "none");
pred_plot <- pred_plot + scale_x_continuous(breaks = c(0,30,60,90));
pred_plot;

ggsave(filename = "figures/section3/prediction.30obs.pdf", plot = pred_plot,
device = "pdf", width = 12, height = 4);

## predict number of patients on each day (not observed patients) -------
predict_it <- function(xT, curr_param, days = (length(y) - obs_length)){
  itpred <- rep(NA, days);
  lambda <- get_lambda(curr_param[3], curr_param[4]);
  gamma <- get_lambda(curr_param[5], curr_param[6]);
  alpha_t <- sis_get_alpha_full_cpp(xT,lambda, gamma);
  for(d in 1 : days){
    agent_states <- (runif(N) < alpha_t);
    itpred[d] <- sum(agent_states);
    alpha_t <- sis_get_alpha_full_cpp(agent_states, lambda, gamma);
  }
  return(itpred)
}
