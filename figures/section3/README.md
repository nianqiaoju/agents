### This folder contains scripts to reproduce figures and tables in Section 3.

All the experiments in this section shares the same running simulated data example: 
	the evolution of an infectious disease on 100 heterogenous individuals over the course of 90 days. 
The dataset is simulated by "setup_hetero.R" and saved to "data_setup_hetero.RData."

## Numerical illustrations of the SMC methods.
	* Produce Table 2 with the script compare_variance_logzhat.R 

## Heatmap of marginal likelihood as a function of parameters. 

	We are interested particularly in whether the MLE converges to the data 
	generating parameters as $T$ increases.

	* 'generate_grid_hetero_slopes.R': fix all parameters except 
	$(\beta_{\lambda}^1, \beta_{\gamma}^2)$ at data generating value. 
	Approximate the marginal likelihood for values of $(\beta_{\lambda}^1, \beta_{\gamma}^2)$ on a grid. 
		** 'data_generate_grid_hetero_slopes.RData': saves data.
		** 'plot_likelihood_slopes.R': makes the heatmaps.
	* 'generate_grid_hetero.R': fix all parameters except $(\beta_{\gamma}^1, \beta_{\gamma}^2)$.

## Bayesian inference of the paramters using PMCMC. 
	* setup_mcmc.R: define pmcmc kernel, prior density;
### performs inference to approximate the posterior $p(\theta\mid y_{0:90})$
	* Define the PMCMC algorithms and Gibbs samplers at 
		** run_pmcmc.R, 
		** run_singlesite_gibbs.R, and 
		run_block5_gibbs.R.
	* The posterior samples are saved at 
		** data_run_pmcmc_final.RData, 
		** data_run_pmcmc_randominit.RData, 
		** data_run_singlesite_gibbs_dgpinit.RData, 
		** data_run_singlesite_randominit.RData, 
		** data_run_block5_gibbs_dgpinit.RData, and 
		** data_run_block5_gibbs_randominit.RData.
	* Analyze and visualize the (p)mcmc chains with the script plot_trace_pmcmc.R.
	* Parameter inference on the individual reproductive numbers at posterior_individual_reproductive_number.R.

### posterior predictive of $p(y_{31:90} \mid y_{0:30})$

	* Generate posterior samples of $p(\theta,X \mid y_{0:30})$ by
		** running the script run_pmcmc_shortobs.R, and 
		** setting the approximate length of observations.
		** saving the samples at data_run_pmcmc_prediction_30obs_dgpinit.RData.
	* Generate posterior predictive samples:
		** using the script simulate_posterior_prediction.R and
		** plotting the posterior preditive trajectories.

## Compare Gibbs with PMCMC:
	* Run the Gibbs samplers, initialized from the DGP and run the prior, with scripts:
		** run_block5_gibbs.R, and 
		** run_singlesite_gibbs.R.
	* The posterior sampler are saved at:
		** data_run_block5_gibbs_randominit.RData,
		** data_run_bloick_gibbs_dgpinit.RData,
		** data_run_singlesite_gibbs_randominit.RData,
		** data_run_singlesite_gibbs_dgpinit.RData.ins;
	* compare_gibbs_pmcmc.R generates histograms for some paramters. The plot is included in the appendix.
