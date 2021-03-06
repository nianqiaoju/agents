% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcondbern.R
\name{rcondbern}
\alias{rcondbern}
\title{Conditional Bernoulli distribution}
\usage{
rcondbern(sum_x, alpha, exact = TRUE, num_sample = NULL)
}
\arguments{
\item{sum_x}{sum of the Bernoullis}

\item{alpha}{vector of success probabilities}

\item{exact}{binary, if TRUE, then use id-check sampling, otherwise use the metropolis algorithm. Default is TRUE.}

\item{num_sample}{default is NULL.}
}
\value{
sample, a length N binary vector, with 1 at the chosen entry
}
\description{
Gives one sample from the conditional Bernoulli distribution, either with exact method or with MCMC. 
This function is a wrapper for idchecking_cpp and metropolis_condbern_cpp
}
