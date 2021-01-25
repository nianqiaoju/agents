#' Compute the TV distance between Poisson-Binomial distribution and its translated Poisson approximation
#' @param alpha vector of probabilities for the Poisson-Binomial distribution
#' @return total variation distance
#' @export
compute_tv_poisbinom_translatedpois <- function(alpha){
  n <- length(alpha)
  p_exact <- rep(0, n+2)
  p_tp <- rep(0, n+2)
  mu <- sum(alpha) 
  sigma2 <- mu - sum(alpha**2)
  translation <- floor(mu - sigma2)
  rate <- mu - translation 
  p_exact[1:(n+1)] <- exp(logdpoisbinom(alpha))
  p_tp[1:(n+1)] <- sapply(0:n, function(x) dpois(x = x - translation, lambda = rate, log = F))
  p_tp[n+2] <- 1 - ppois(q = n - translation, lambda = rate, log = F) ## this term deals with the unequal support
  tv <- sum(abs(p_tp - p_exact)) / 2
  return(tv)
}
