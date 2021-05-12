#' @title Metropolis algorithm targeting the conditional Bernoulli distribution 
#' @description The MCMC transition kernel for the conditional Bernoulli distribution. 
#' Given the current state x, 
#' it proposes to switch a randomly chosen 1 with a randomly chosen 0. 
#' @param sum_x sum of the samples 
#' @param alpha vector of probabilities, length N
#' @param config list of configurations, save_history, num_sample
#' @return either xx (a vector of size N) or  xx_matrix (a binary matrix of size num_samples by N)
#' @export


metropolis_condbern <- function(sum_x, alpha,num_sample = 10, save_history = FALSE){
  N_ <- length(alpha)
  ## initalize entries with 1 and entries with 0  
  n1 <- sample.int(n = N_, size = sum_x) ## there is a 1 at these locations
  n0 <- setdiff(c(1:N_), n1) ## there is a 0 at these locations
  if (save_history){
    xx_matrix <- matrix(0, nrow = num_sample, ncol = N_)
  }
  for (iter in 1:num_sample){
      j1 <- ceiling(runif(1) * sum_x) ## index in the set n1
      j0 <- ceiling(runif(1) * (N_ - sum_x)) ## index in the set n0
      ## xx[n1[j1]] = 1 and xx[n0[j0]] = 0  
      accept_probability <-  alpha[n0[j0]] * ( 1 - alpha[n1[j1]]) / alpha[n1[j1]] / ( 1 - alpha[n0[j0]])
      if (runif(1) < accept_probability){
        ## swap the locations 
        temp <- n1[j1]
        n1[j1] <- n0[j0]
        n0[j0] <- temp
      }
      if (save_history){
        xx_matrix[iter, n1] <- 1
      }
  }
  if (save_history){
    return(xx_matrix)
  } else {
    xx <- rep(0, N_)
    xx[n1] <- 1
    return(xx)
  }
}
