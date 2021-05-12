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
# 
#   ## compute all the q(h,n) from n = N to n = 1 and h = 0,1,..., N 
#   N <- length(alpha)
#   previous_log_densities <- rep(-Inf, N + 1)
#   current_log_densities <- rep(-Inf, N + 1)
#   
#   # initialize the recursion, start with q(0,N) and q(1,N)
#   previous_log_densities[0 + 1] <- log(1 - alpha[N])
#   previous_log_densities[1 + 1] <- log(alpha[N])
#   if (N == 1){## if alpha is a scalar
#     return(previous_log_densities)
#   }
#   for (n_complement in 1:(N-1)){
#     ## vector under consideration is alpha[n : N] which has length n_complement
#     n <- N - n_complement
#     current_log_densities[0 + 1] <- previous_log_densities[0 + 1] + log(1 - alpha[n])
#     for (h in 1 : (n_complement + 1)){ ## the partial sum cannot exceed length of the vector currently under inspection
#       ls1 <- log(alpha[n]) + previous_log_densities[h - 1 + 1]
#       ls2 <- log(1 - alpha[n]) + previous_log_densities[h + 1]
#       maxls <- max(ls1,ls2)
#       current_log_densities[h + 1] <- log(exp(ls1 - maxls) + exp(ls2 - maxls)) + maxls 
#     }
#     previous_log_densities <- current_log_densities
#     current_log_densities <- rep(-Inf, N + 1)
#   }
#   return(previous_log_densities)
# }

## Note: all the log(x) and log(1-x) could be computed once and for all (i.e. N times instead of N^2)

# logpoisbinom <- function(s, prob){
#   pmf <- dpoisbinom(s, prob)
#   negative <- (pmf < 0) 
#   pmf[negative] <- 0
#   logpmf <- log(pmf)
#   return(logpmf)
# }
