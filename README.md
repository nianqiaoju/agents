This repository contains code for the paper `Sequential Monte Carlo algorithms for agent-based models of disease transmission' by Nianqiao (Phyllis) Ju, Jeremy Heng and Pierre Jacob. 

## Installation
The package can be installed from R via:
``` r
# install.packages("devtools")
devtools::install_github("nianqiaoju/agents")
```
It depends on the packages Rcpp, RcppEigen, lubridate, which can be
installed via:

``` r
install.packages(c("Rcpp", "RcppEigen"))
```

Additionally you might want to install other packages, to help with
parallel computation and to compare run-time of algorithms:

``` r
install.packages(c("doParallel", "doRNG", "tictoc"))
```

and to help with manipulating results and plotting:

``` r
install.packages(c("dplyr", "ggplot2", "ggtheme"))
```

although these packages are not strictly required.