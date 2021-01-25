### scripts to reproduce figures & tables of Section 2

This folder contains scripts to reproduce the numerical
results of Section 2 ("Static agent-based model"), regarding
the translated Poisson approximation to the Poisson Binomial distribution,
the MCMC sampler for Conditional Bernoulli distributions, and 
a numerical comparison of various MCMC samplers
to estimate the parameters and latent variables.

### Instructions

The packages {"agents", "dplyr", "gridExtra", "tictoc", "reshape"} must be installed.

Then the following scripts must be executer, in the specified order,
from the working directory set to the root folder of the 'agents' package.

Table 1 
## set up features and data for the static model
* setup_random_alpha.R
## define random walk MH and Gibbs sampling functions and tuning parameters
* tune_mh.R
## approximate the Monte Carlo error, by running a long MH-exact chain and independent repeats of MH-exact, MH-tp, PMMH and Gibbs
* sample_posterior_long_n1000.R
* repeat_posterior_n1000.R
* compute_mse.R

______
the following files have been moved to backup, 
since the figures have been removed from the paper.

Figure 1 a [Note: now uses the C++ implementation of logdpoisbinom]
#### this we can potentially remove
* plot_transpoi_tvdistance.R

Figure 1 b [Note: now uses the C++ implementation of the coupled MCMC scheme for CBernoulli]
#### this we can remove as it is relegated to the other paper
* plot_condber_mixingtime.R

