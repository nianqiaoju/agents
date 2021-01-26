### This folder (/figures/section3) contains scripts to reproduce figures and tables in section 3.

All the experiments in this section shares the same running simulated data example: the evolution of an infectious disease on 100 heterogenous individuals over the course of 90 days. 
"setup_hetero.R" generates the simulated dataset and saves it at "data_setup_hetero.RData."

Experiments in this section serve three purposes.
## Properties of the three SMC algorithms, in terms of run-time, effective sample size and variance of marginal likelihood estimators. 
	* compare_variance_logzhat produces Table 2. 
## Heatmap of marginal likelihood as a function of parameters. We are interested particularly in whether the MLE converges to the data generating parameters as $T$ increases.
	* 'generate_grid_hetero_slopes.R': fix all parameters except $(\beta_{\lambda}^1, \beta_{\gamma}^2)$ at data generating value. Approximate the marginal likelihood for values of $(\beta_{\lambda}^1, \beta_{\gamma}^2)$ on a grid. 
		** 'data_generate_grid_hetero_slopes.RData': saves data.
		** 'plot_likelihood_slopes.R': makes the heatmaps.
	* 'generate_grid_hetero.R': fix all parameters except $(\beta_{gamma}^1, \beta_{\gamma}^2)$.

## Bayesian inference of the paramters using PMCMC. 
	* setup_mcmc.R: define pmcmc kernel, prior density;
### performs inference to approximate the posterior $p(\theta\mid y_{0:90})$
	* run_pmcmc.R, run_singlesite_gibbs.R and run_block5_gibbs.R defines the PMCMC algorithm and the Gibbs samplers.  
	* The posterior samples are saved at data_run_pmcmc_final.RData, data_run_pmcmc_randominit.RData, data_run_singlesite_gibbs_dgpinit.RData, data_run_singlesite_randominit.RData, data_run_block5_gibbs_dgpinit.RData, and data_run_block5_gibbs_randominit.RData.
	* Analyze and visualize the (p)mcmc chains with the script plot_trace_pmcmc.R.
	* Parameter inference on the individual reproductive numbers at posterior_individual_reproductive_number.R.

### posterior predictive of $p(y_{31:90} \mid y_{0:30})$'run_pmcmc_shortobs.R': ;
	* generate posterior samples of $p(\theta,X \mid y_{0:30})$ by running the script run_pmcmc_shortobs.R and setting the approximate length of observations.
	* generate posterior predictive samples with the script run_posterior_prediction.R and save plots.

## Compare Gibbs with PMCMC:
	* run_block5_gibbs.R
	* run_singlesite_gibbs.R
	* data_run_block5_randinit.RData and data_run_singlesite_gibbs_randinit.RData contains the resulting chains;
	* compare_gibbs_pmcmc.R generates histograms for some paramters. The plot is included in the appendix.
