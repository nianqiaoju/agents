% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sis_apf.R
\name{sis_apf}
\alias{sis_apf}
\title{Fully adapted auxiliary Particle Filter for SIS model with Population Observation}
\usage{
sis_apf(y, model_config, particle_config)
}
\arguments{
\item{y}{a vector of length (T+1)}

\item{model_config}{a list containing:
\itemize{
\item 'N': size of the population
\item 'alpha0' : initial infection probability
\item 'lambda': infection rate 
\item 'gamma': recovery rate
\item 'adjacency_matrix_b': network structure
\item 'rho': reporting rate
\item }}

\item{particle_config}{a list containing:
\itemize{
\item 'num_particles' : number of particles for the bootstrap particle filter
\item 'ess_threshold' : if effective sample size drops below the threshold, then perform a resample step. ess_threshold = 1 means resampling at every step.
\item 'save_particles': binary
\item 'clock' : binary, default to FALSE. If clock = TRUE, then we will use a stopwatch to document its Sys.time()
\item 'save_genealogy': binary
\item 'verbose':
\item 'exact': binary, if TRUE, then uses exact PoisBinom and exact CondBern, if FALSE then use translated Poisson and MCMC for CondBern. MUST HAVE!
\item 'num_mcmc": required if exact = FALSE;
}}
}
\value{
A list containing 
\itemize{
\item 'log_incremental_likelihood' :  
\item 'log_marginal_likelihood' :
\item 'log_final_likelihood' :
\item 'particles': [N, (T+1), num_particles] array storing the particles: 
\item 'eves' : the unique ancesters of each particle in the genealogy structure, if save_genealogy = TRUE;
\item 'ess' : effective sample size at each step
\item 'runtime' : elapse time at each step, if clock = TRUE;
\item 'totaltime' : elapse time of running the whole particle filter, if clock = TRUE;
}
}
\description{
Fully adapted auxiliary Particle Filter for SIS model with Population Observation
}
