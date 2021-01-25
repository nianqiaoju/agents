## this script analyses the posterior samples
## it studies the individual reproductive numbers, R_0^n = \lambda^n / \gamma^n;

rm(list = ls());
library(dplyr);
library(ggplot2);
source("inst/set_plots.R");

## load data -----
# load("figures/section3/data_run_pmcmc_final.RData");
load("figures/section3/data_run_pmcmc_randominitial.RData");
pmcmc_accept_chain_dgpinit <- accept_chain;
pmcmc_lpost_chain_dgpinit <- lpost_chain;
pmcmc_param_chain_dgpinit <- chain;

dim(pmcmc_param_chain_dgpinit)

nburn <- burn_in;
nmcmc <- num_mcmc;

## functions to credible sets and other summary stats ----
get_confset <- function(vec){
  quantile(vec, probs = c(0.025, 0.5,0.975));
}

get_mean_sd <- function(vec){
  vec <- vec[!is.na(vec)];
  c(mean(vec) - 2 *  sd(vec), mean(vec), mean(vec) + 2 * sd(vec));
}


## print marginal summaries -----
for(id in 1 : 7){
  cat("dgp is ", dgp_parameters[id],"\n");
  cat(get_confset(pmcmc_param_chain_dgpinit[nburn : nmcmc, id]), "\n");
}

## helper functions for individual reproductive numbers
get_lambda <- function(b1,b2){
  agents::expit(b1 * features[1,] + b2 * features[2,]);
}

get_r0_all <- function(parameters){
  ## lambda / gamma
 get_lambda(parameters[3], parameters[4]) / get_lambda(parameters[5], parameters[6]);
}

get_r0_agent <- function(parameters, agent_id){
  agents::expit(parameters[3] * features[1,agent_id] + parameters[4]* features[2, agent_id]) / agents::expit(parameters[5] * features[1, agent_id] + parameters[6] * features[2, agent_id]);
}

## individual reproductive numbers ----
rn_df <- data.frame(id = integer(N), lbd = double(N), med = double(N), ubd = double(N), dgp = double(N));
rn_df$id <- 1 : N;
rn_df$dgp <- log(get_r0_agent(dgp_parameters)); 
for(id in 1 : N){
  rn_df[id, c(2:4)] <- get_confset(apply(pmcmc_param_chain_dgpinit[nburn : nmcmc, ], 1, function(param) log(get_r0_agent(param, id))));
}
head(rn_df);
rn_df <- rn_df %>% arrange(dgp);
head(rn_df);

## creat a plot for the posterior distribution

rn_plot <- ggplot(rn_df, aes(x = 1 : N, y = med)) ;
rn_plot <- rn_plot + geom_errorbar(aes(ymin=lbd, ymax=ubd), width= 1.5 , linetype = 1,  color = cbPalette[3]);
rn_plot <- rn_plot + geom_point(aes(y = dgp), color = "black") + geom_point( color = "blue");
rn_plot <- rn_plot + labs(x = "sorted agent index", y = expression(log(R[0]^n))) + theme(legend.position = "none"); 
rn_plot

# ## save g to pdf 
ggsave(filename = "figures/section3/individual.reproductive.number.pdf", plot = rn_plot,
       device = "pdf", width = 6, height = 6);

## the distribution of reproductive numbers -----

## compute logr0 for each agent at the dgp
get_log_r0 <- function(parameters) log(get_r0_all(parameters));
lr0_dgp <- get_log_r0(dgp_parameters);


lpt <- seq(-10,9, 0.75);
rpt <- c(lpt[-1], 9 + .75);
mpt <- (lpt + rpt) * 0.5;
nbins <- length(lpt);
df <- data.frame(lpt = lpt, mpt = mpt, rpt = rpt)
df$upper_envelope <- df$lower_envelope <- df$median_evelope <- df$dgp_freq <- double(nbins)


## compute the posterior samples of logr0
## and place them into bins
posterior_lr0 <- apply(pmcmc_param_chain_dgpinit[nburn : nmcmc, ], 1, get_log_r0);
for(ibin in 1 : dim(df)[1]){
  ## get the percerntage of points in each bin from data generating values
  df$dgp_freq[ibin] <- mean(lr0_dgp > df$lpt[ibin] & lr0_dgp < df$rpt[ibin]);    
  ## get the percentage of points in each bin from posterior parameters
  posterior_lr0_in_bin <- apply(posterior_lr0, 2, function(x) mean(x > df$lpt[ibin] & x < df$rpt[ibin])) ;
  df$lower_envelope[ibin] <- quantile(posterior_lr0_in_bin, 0.025);
  df$upper_envelope[ibin] <- quantile(posterior_lr0_in_bin, 1 - 0.025);
  df$median_evelope[ibin] <- median(posterior_lr0_in_bin);
}
head(df)

## organize df into df_hist for ggplot
df_hist <- data.frame(x = double(2 * (nbins - 1)), y = double(2 * (nbins - 1)), z = double(2 *  (nbins - 1)), lbd = double(2 *  (nbins - 1)), ubd = double(2 *  (nbins - 1)) );
for(ibin in 1 : (nbins - 1)){
  df_hist[2 * ibin - 1, ] <-  c(df$lpt[ibin], df$dgp_freq[ibin], df$median_evelope[ibin],  df$lower_envelope[ibin], df$upper_envelope[ibin]) ## left end
  df_hist[2 * ibin, ] <- c(df$lpt[ibin + 1], df$dgp_freq[ibin], df$median_evelope[ibin], df$lower_envelope[ibin], df$upper_envelope[ibin])
}

## plot the histograms
hist_plot <- ggplot(df_hist, aes(x = x, y = y, ymin = lbd, ymax = ubd));
hist_plot <- hist_plot +  geom_ribbon(alpha = 0.2, fill = cbPalette[3]) + 
  geom_line(aes(y = lbd), color = cbPalette[3], linetype = 2, lwd = 1) + 
  geom_line(aes(y = ubd), color = cbPalette[3], linetype = 2, lwd = 1) +
  geom_line(lwd = 1, color = "black") +
  geom_line(aes(y = z), color = "blue", lwd = 1) ;
hist_plot <- hist_plot + labs(x = expression(log(R[0]^n)), y = "proportion");
# hist_plot <- hist_plot + theme(axis.text.y = element_text(angle = 90));
hist_plot

ggsave(filename = "figures/section3/profile.reproductive.number.pdf", plot = hist_plot,
       device = "pdf", width = 6, height = 6);

## plot posterior histogram envelops of the individual r0 -------
### 95% credible envelopes for log(lambda / gamma)
plot(df_hist$x, df_hist$y, xlim = range(lpt), ylim = c(0,1) , type = "l",
     main = "95% credible envelopes for log(lambda/ gamma)", 
     xlab = "log(lambda / gamma)",
     ylab = "frequency")
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(df_hist$x,rev(df_hist$x)),c(df_hist$lbd,rev(df_hist$ubd)),col = c1, border = FALSE)
lines(df_hist$x, df_hist$y, lwd = 2)
lines(df_hist$x, df_hist$z, col = "red", lwd = 2)
#add red lines on borders of polygon
lines(df_hist$x, df_hist$lbd, col="red",lty=2)
lines(df_hist$x, df_hist$ubd, col="red",lty=2)

