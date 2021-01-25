rm(list = ls())
load("figures/section3/data_compare_ess_3regimes.RData")
library(ggplot2)
library(ggthemes)
agents::setmytheme()

# change to percentage scale
ess_df$ess <- ess_df$ess * 100

# change label of regimes
ess_df$regime <- factor(ess_df$regime, levels = c(1,2,3), labels = c("original", "small", "large"))

# plot
g <- ggplot(ess_df, aes(x = step, y = ess, colour = method)) + geom_line() + geom_point()
g <- g + scale_color_colorblind()
g <- g + facet_grid(rows = vars(regime)) 
g <- g + xlab("time") + ylab("ESS%") + xlim(c(0, 90))
g

# ## save g to pdf 
ggsave(filename = "~/Dropbox/AgentBasedModels/paper/sis.compare.ess.3regimes.pdf", plot = g,
       device = "pdf", width = 12, height = 6);
