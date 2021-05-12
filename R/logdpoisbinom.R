#' @title Poisson Binomial distribution
#' @description Density function for Poisson binomial distribution with parameters alpha, 
#' which is the distribution of the sum of Bernoulli random variables, each with a different
#' probability \eqn{P(X_i = 1) = \alpha_i}.
#' @param alpha vector of N probabilities (numbers between 0 and 1).
#' @return the entire density function in a vector of length N+1.
#' @export

logdpoisbinom <- function(alpha){
  logdpoisbinom_cpp(alpha)
}

