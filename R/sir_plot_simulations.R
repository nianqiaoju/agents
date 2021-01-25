#' @export
plot_sir_simulations <- function(result,t){
  g <- result$graph
  layout <- result$layout
  hidden_status<- result$hidden_status
  #y <- result$observations
  colors = as.character(factor(hidden_status[,t],
                               levels = c('S','I','R'), 
                               labels = c("chartreuse3","indianred1","gold")))
  plot(g, layout = layout,
       edge.arrow.size = .5,  
       vertex.color = adjustcolor(colors , alpha.f = 0.8), 
       vertex.label = NA, vertex.size = 5)   #, main = paste(sim[[4]],"Network, t =", t)) 
  legend(x = 1, y = 0,
         inset=c(0,0),
         c('Susceptible','Infected','Recovered'),
         pch = 21, 
         pt.bg = c("chartreuse3","indianred1","gold"), 
         bty = "n", 
         pt.cex = 1.5)
}
