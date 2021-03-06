% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sis_forward_backward.R
\name{sis_forward_backward}
\alias{sis_forward_backward}
\title{Forward backward algorithm for the SIS model}
\usage{
sis_forward_backward(y, model_config, num_samples = 1, all_xx = NULL)
}
\arguments{
\item{y}{population observations}

\item{model_config}{a list containing
\itemize{
\item 'N' : population size
\item 'alpha0'
\item 'lambda'
\item 'gamma'
\item 'rho'
\item 'adjacency_matrix_b': network structure
}}

\item{num_samples}{number of backward samples needed, set to 1 if not specifed. If num.samples = 0, then no backward sample is computed.}

\item{all_xx}{a 2**N by N matrix, representing the sample space, pre-calculated to save time}
}
\value{
a list containing 
\itemize{
\item 'log marginal likelihood and one backward sample' 
\item 'bs' : backward samples in a [N * num_observations * num_samples] array
}
}
\description{
Forward backward algorithm for the SIS model
}
