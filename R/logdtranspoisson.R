#' @title translated Poisson distribution
#' @description density approximation of the Poission binomial distribution in log scale
#' @param xx vector of integers
#' @param mu mean 
#' @param sigmasq variance
#' @return log densities
#'
#' @export

logdtranspoisson <- function(xx, mu, sigmasq){
  gammaa <- (mu - sigmasq) %% 1
  translation <- floor(mu - sigmasq)
  if (min(xx) >= translation){
    logpmf <- dpois(x = xx - translation, lambda = sigmasq + gammaa, log = TRUE)
  } else {
    logpmf <- rep(-Inf, length(xx))
    index <- (xx >= translation)
    logpmf[index] <- dpois(x = xx[index] - translation, lambda = sigmasq + gammaa, log = TRUE)
  }
  return(logpmf)
}