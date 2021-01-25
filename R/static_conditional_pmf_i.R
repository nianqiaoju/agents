#' @title Conditional PMF of I in the static model 
#' @description computes the marginal likelihood observing y given alpha in the static case and computes the conditional pmf of each i 
#' @param alpha a length N vector of probabilities
#' @param y observation (0 <= y <= N)
#' @param rho reporting rate
#' @param N length of alpha
#' @return a list containing 
#' \itemize{
#' \item cpmf:conditional pmf of i given y as a vector of probabilities
#' \item llik:log-likelihood of y 
#' }
#' @export
#' 
#' 
static_conditional_pmf_i <- function(y, alpha, rho, N){
  cpmf <- rep(NA, N - y  + 1)
  if (y==N){
    cpmf[1] <- 1
    llik <- sum(log(alpha)) + N * log(rho)
    return(list(cpmf = cpmf, llik = llik))
  }else{
    log_cpmf <- logdpoisbinom(alpha) + dbinom(x = y, size = 0 : N, prob = rho, log = TRUE)
    cpmf <- lw.normalize(log_cpmf)
    llik <- lw.logsum(log_cpmf) ## should be the same as logdpoisbinom(alpha * rho)[y + 1]
    return(list(cpmf = cpmf[(y + 1) : (N + 1)], llik = llik))
  }
}
