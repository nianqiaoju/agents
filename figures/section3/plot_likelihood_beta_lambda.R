rm(list = ls());
### no need to set.seed since there is nothing random

load("figures/section3/data_generate_grid_hetero.RData");

library("ggplot2");
library(dplyr);
library("scales")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##
cat("[reminder: data generating parameters are: bl1 = ", dgp_bl1,",bl2 = ",dgp_bl2,"]\n");

# llik <- subset(llik, llik != 0); ## DEBUG CODE


## mutate the likelihood a bit for better contour plots
## find mle for each value of t
# llik <- llik - maxllik for each t, so that mle values are contour plot is easier to see


mle_summary <- data.frame(llik_df %>% group_by(t) %>% summarise(maxllik = max(llik), bl1 = bl1[which(llik == maxllik)], bl2 = bl2[which(llik == maxllik)]));
for(step in time_grid){
  llik_df$llik[which(llik_df$t == step)] <- llik_df$llik[which(llik_df$t == step)] - mle_summary$maxllik[which(mle_summary$t == step)];
}
## change t to factor
llik_df$t <- factor(llik_df$t, levels = time_grid, labels = paste("T =", time_grid));
mle_summary$t <- factor(mle_summary$t, levels = time_grid, labels = paste("T =", time_grid));

## generate contour plot
setmytheme <- function(){
  library(ggplot2)
  library(ggthemes)
  theme_set(theme_bw())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
               axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               legend.position = "bottom")
}
setmytheme();
g <- ggplot(data = llik_df, aes(x = bl1, y = bl2));
g <- g + geom_contour(aes(z = llik, colour = stat(level)), bins = 10);
g <- g + geom_point(aes(x = dgp_bl1, y = dgp_bl2), colour = "red", size = 4);
g <- g + facet_grid(. ~ t);
g <- g + geom_point(data = mle_summary, aes(x = bl1, y = bl2, group = t), colour = "green", size = 4);
g <- g + scale_colour_distiller(palette = "YlGn", direction = 1) + guides(colour = guide_colorbar(title = "log-likelihood", barwidth = 20));
g;

llik_df$llik[which(llik_df$llik < - 100)] <- -100

heatmap <- ggplot(data = llik_df, aes(x = bl1, y = bl2, group = t));
heatmap <- heatmap + facet_grid( . ~ t)+ geom_tile(aes(fill = llik)) ;
heatmap <- heatmap + geom_point(aes(x = dgp_bl1, y = dgp_bl2), colour = "black", size = 4) + 
  geom_point(data = mle_summary, aes(x = bl1, y = bl2, group = t), colour = "red", size = 4) ;
heatmap <- heatmap + scale_fill_gradientn(colours=c("white", cbPalette[8], cbPalette[3], "darkgreen", cbPalette[5]),
                                          values=rescale(c(-100, -20, -10, -5,0)),
                                          breaks=c(-100, -20, -10, -5,0),
                                          labels=c(-100, -20, -10, -5,0),
                                          guide="colorbar")
heatmap <- heatmap + guides(fill = guide_colorbar(title = "log-likelihood", barwidth = 40), ticks =  TRUE, labels = TRUE) ;
heatmap <- heatmap + xlab(expression(beta[lambda]^1)) + ylab(expression(beta[lambda]^2))
heatmap

ggsave(filename = "~/Dropbox/AgentBasedModels/paper/sis.likelihood.heatmap.hetero.pdf", plot = heatmap,
       device = "pdf", width = 12, height = 6);
