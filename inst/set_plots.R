## ggplot2 set up 
## color blind pallete for ggplot
cbPalette <- c("#101010", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## set theme to blackwhite and determine margins
agents::setmytheme();


## some colors for overlapping histograms
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue");
c3 <- rgb(173,216,230,max = 255, alpha = 100, names = "lt.blue");
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink");

library(ggthemes);
