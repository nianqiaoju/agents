rm(list = ls());
library(ggplot2);
library(agents);
source("inst/set_plots.R")
load("figures/section3/data_run_pmcmc_prediction_30obs_dgpinit.RData");
set.seed(2020);

## predict number of observed patients on each day ----
sis_predict <- function(xT, curr_param, days = (length(y) - 1 - obs_length)){
  ypred <- rep(NA, days);
  lambda <- get_lambda(curr_param[3], curr_param[4]);
  gamma <- get_lambda(curr_param[5], curr_param[6]);
  alpha_t <- sis_get_alpha_full_cpp(xT,lambda, gamma);
  for(d in 1 : days){
    agent_states <- (runif(N) < alpha_t);
    ypred[d] <- rbinom(1, sum(agent_states), curr_param[7]);
    alpha_t <- sis_get_alpha_full_cpp(agent_states, lambda, gamma);
  }
  return(ypred)
}

nmcmc <- 105000
nburn <- 55000

ypred_samples <- matrix(NA, (length(y) - obs_length - 1), nmcmc -nburn);
for(isample in (nburn + 1) : nmcmc){
  ## simulate a trajectory
  ypred_samples[, isample - nburn] <- sis_predict(xT_chain[isample], param_chain[isample,]);
}

## how many of the predictions leads to no breakout?
no_outbreaks <- function(ypred){
  sum(ypred == 0) >= 1;
}

## ~43 of the predicted trajectories are not outbreaks
mean(apply(ypred_samples, 2 , no_outbreaks));


sample_id <- which(ypred_samples[1,] < 20 & ypred_samples[1,] > 5 & ypred_samples[10,] < 5);
thining_id <- c(1013, 940, sample_id, seq(1, dim(ypred_samples)[2], length.out = 30))
## plot data and predictions
plot( x= 0 : (length(y) - 1) , y = y, main = "observations and predictions", ylim = c(0,80), xlim = c(0, length(y) - 1),
      pch = 1, 
      xlab = "time",
      ylab = "observations")
## plot trajectories after thinning
for(isample in thining_id){
  lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_samples[,  isample], col = "grey75", lwd = 1);
}
# for(isample in sample_id){
#   lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_samples[,  isample], col = "grey75", lwd = 1);
# }
lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = y[(obs_length  + 2 ) : length(y)], col = "black", lwd = 2)
points( (obs_length  + 1 ) : (length(y) - 1 ), y = y[(obs_length  + 2 ) : length(y)], col = "black",
        pch = 1);
points( 0 : obs_length, y = y[1 : (obs_length + 1)], pch = 19);

## plot data and confidence curves
## get median and quantiles y for each day
get_quantiles <- function(xx){
  quantile(xx, probs = c(0.025,0.5, 1 - 0.025));
}

y_pred_quantiles <- apply(ypred_samples, 1, get_quantiles);

lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = y_pred_quantiles[1, ], col = "red", lty = 2, lwd = 2);
lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = y_pred_quantiles[2, ], col = "red", lty = 2, lwd = 2);
lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = y_pred_quantiles[3, ], col = "red", lty = 2, lwd = 2);

## plot confidence band for those trajactories that are outbreaks
# ypred_samples_outbreaks <- ypred_samples[, !apply(ypred_samples, 2 , no_outbreaks) ]
# ypred_outbreaks_quantiles <- apply(ypred_samples_outbreaks, 1, get_quantiles);
# 
# lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_outbreaks_quantiles[1, ], col = "yellow", lwd = 2, lty = 2);
# lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_outbreaks_quantiles[2, ], col = "yellow", lwd = 2, lty = 2);
# lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_outbreaks_quantiles[3, ], col = "yellow", lwd = 2, lty = 2);

get_nozero_quantile <- function(xx){
  quantile(xx[!xx==0], probs = c(0.025,0.5,1-0.025));
}

ypred_nozero_quantiles <- apply(ypred_samples, 1, get_nozero_quantile);

lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_nozero_quantiles[1, ], col = "blue", lwd = 2, lty = 2);
lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_nozero_quantiles[2, ], col = "blue", lwd = 2, lty = 2);
lines(x = (obs_length  + 1 ) : (length(y) - 1 ), y = ypred_nozero_quantiles[3, ], col = "blue", lwd = 2, lty = 2);


## a dataframe for the trajectories
traj_df <- data.frame(x = integer(length(y) - obs_length - 1), y = integer(length(y) - obs_length - 1));
traj_df$x <-  (obs_length  + 1 ) : (length(y) - 1 )
head(traj_df);
## a dataframe for the quantiles
quantile_df <- data.frame(x = (obs_length  + 1 ) : (length(y) - 1 ), 
                          ypredmin = y_pred_quantiles[1, ], 
                          ypredmed = y_pred_quantiles[2, ],
                          ypredmax = y_pred_quantiles[3, ],
                          yprednozeromin = ypred_nozero_quantiles[1, ],
                          yprednozeromed = ypred_nozero_quantiles[2, ], 
                          yprednozeromax = ypred_nozero_quantiles[3, ]);

data_df <- data.frame(x = 0 : (length(y) - 1 ), y = y, observed = c(rep(TRUE, 31), rep(FALSE, 60)));
data_df$observed <- factor(data_df$observed);

pred_plot <- ggplot(data_df, aes(x = x, y = y));
pred_plot <- pred_plot + geom_point(size = 3, shape = 21, aes(fill = observed));
for(isample in thining_id){
  traj_df$y <- ypred_samples[,  isample];
  pred_plot <- pred_plot + geom_line(data = traj_df, aes(x = x, y = y), color = "grey80");
}
pred_plot <- pred_plot + geom_line(data = quantile_df, aes(x = x, y = ypredmed), linetype = 2, size = 2, color = "blue") + 
  geom_line(data = quantile_df, aes(x = x, y = ypredmin), linetype = 2, size = 2, color = cbPalette[3]) +
  geom_line(data = quantile_df, aes(x = x, y = ypredmax), linetype = 2, size = 2, color = cbPalette[3]);
# pred_plot <- pred_plot + geom_line(data = quantile_df, aes(x = x, y = yprednozeromed), linetype = 2, size = 2, color = cbPalette[2]) + 
#   geom_line(data = quantile_df, aes(x = x, y = yprednozeromin), linetype = 2, size = 2, color = cbPalette[2]) +
#   geom_line(data = quantile_df, aes(x = x, y = yprednozeromax), linetype = 2, size = 2, color = cbPalette[2]);
pred_plot <- pred_plot + geom_line(data = data_df[(obs_length + 2) : length(y), ]);
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
