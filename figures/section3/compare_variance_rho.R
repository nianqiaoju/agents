rm(list = ls());
library(dplyr)
library(ggplot2)
library(ggthemes)
setmytheme()
load("figures/section3/data_generate_grid_rho_fix_beta.RData");

grid_rho_df$method <- factor(rep(c("cSMC1", "cSMC2", "BPF", "APF"), 900), levels = c("cSMC1", "cSMC2", "BPF", "APF")) 

grid_rho_summary <- grid_rho_df %>% group_by(rho, method) %>% summarise(varlogz = var(llik));
tail(grid_rho_summary)
class(grid_rho_summary)

grid_rho_summary$method <- factor(grid_rho_summary$method, levels = c("BPF", "APF", "cSMC1", "cSMC2"))
g <- ggplot(grid_rho_summary, aes(x = rho, y = varlogz, group = method, colour = method)) + geom_line()
g <- g + scale_color_colorblind() + xlab(expression(rho)) + ylab("Variance") 
g <- g + scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1))
g <- g + scale_y_log10()
g

# save plot
ggsave(filename = "~/Dropbox/AgentBasedModels/paper/sis.variance.rho.pdf", plot = g,
       device = "pdf", width = 8, height = 6);

# compute efficiency relative to BPF
efficiency_df <- data.frame()
relative_cost <- 1 / 19 
for (i in 1:length(rho_grid)){
  rho <- rho_grid[i] 
  irow <- 4 * (i-1)
  csmc1_result <- grid_rho_summary[irow + 1,]
  csmc2_result <- grid_rho_summary[irow + 2,]
  bpf_result <- grid_rho_summary[irow + 3,]
  apf_result <- grid_rho_summary[irow + 4,]

  # csmc1
  efficiency_df <- rbind(efficiency_df, data.frame(rho = rho, 
                                                   method = csmc1_result$method, 
                                                   efficiency = relative_cost * bpf_result$varlogz / csmc1_result$varlogz))  

  # csmc2
  efficiency_df <- rbind(efficiency_df, data.frame(rho = rho, 
                                                   method = csmc2_result$method, 
                                                   efficiency = relative_cost * bpf_result$varlogz / csmc2_result$varlogz))  
  
  # apf
  efficiency_df <- rbind(efficiency_df, data.frame(rho = rho, 
                                                   method = apf_result$method, 
                                                   efficiency = relative_cost * bpf_result$varlogz / apf_result$varlogz))  
  # bpf
  efficiency_df <- rbind(efficiency_df, data.frame(rho = rho, 
                                                   method = bpf_result$method, 
                                                   efficiency = 1))  
  
}

# plot efficiency relative to BPF
efficiency_df$method <- factor(efficiency_df$method, levels = c("BPF", "APF", "cSMC1", "cSMC2"))
g <- ggplot(efficiency_df, aes(x = rho, y = efficiency, group = method, colour = method)) + geom_line()
g <- g + scale_color_colorblind() + xlab(expression(rho)) + ylab("Relative efficiency") 
g <- g + scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1))
# g <- g + scale_y_log10()
g

# save plot
ggsave(filename = "~/Dropbox/AgentBasedModels/paper/sis.efficiency.rho.pdf", plot = g,
       device = "pdf", width = 8, height = 6);