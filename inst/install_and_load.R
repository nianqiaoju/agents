# run this file to install all the required packages to run the package agent

list.of.packages <- c("ggplot2","dplyr","gridExtra","tictoc", "ggthemes")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


## load the packages
library(ggplot2)
library(agents)
library(gridExtra)
library(dplyr)
library(ggthemes)
# library(shiny)
# library(igraph)
# library(deSolve)
# library(tictoc)



