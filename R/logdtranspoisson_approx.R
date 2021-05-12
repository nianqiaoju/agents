#' @title translated Poisson approximation 
#' @description density approximation of the Poission binomial distribution in log scale. \\
#' NOTE: the densities do not sum to 1.
#' @param alpha vector of Bernoulli parameters
#' @param xx  default is NULL, when x == NULL, output is the entire density function 
#' @return approximate log densities. 
#'
#' @export

logdtranspoisson_approx <- function(alpha, y = NULL){
  mu <- sum(alpha);
  sigma2 <- mu - sum(alpha**2);
  translation <- floor(mu - sigma2);
  rate <- mu - translation;
  if (is.null(y)) y <- c(0 : length(alpha));
  return(dpois(x =  y - translation, lambda = rate, log = TRUE));
}