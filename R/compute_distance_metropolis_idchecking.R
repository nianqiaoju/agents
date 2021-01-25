#' @title The Metropolis sampler for Conditional Bernoulli distribution
#' @description compute an upper bound for the 2-Wasserstein distance 
#' between Metropolis sampler in metropolis_condbern.R and the Conditional-Bernoulli distribution 
#' @param alpha vector of probabilities for the Poisson-Binomial distribution 
#' @param sum_x sum of all entries
#' @param num_mcmc number of MCMC iterations to run before coupling 
#' @param num_repeat number of repetitions to get an estimate of the upper bound 
#' @param distance can be 'tv' or 'wasserstein'
#' @return a vector of distance upper bounds for different number of iterations
#' @export
#' 
compute_distance_metropolis_to_condbern <- function(alpha, sum_x, num_mcmc, num_repeat, distance = 'tv'){
  if (distance == 'wasserstein'){
    sse_matrix <- matrix(0, nrow = num_mcmc, ncol = num_repeat)
    for (irep in 1:num_repeat){
      sse_matrix[,irep] <- couple_metropolis_to_condbern(alpha, sum_x, num_mcmc, compute_sse = TRUE, verbose = FALSE)$sse_vector
    }
    wass_upper_bounds <- sqrt(rowMeans(sse_matrix))
    return(wass_upper_bounds)
  }else{
    indicator_of_coupling <- matrix(0, nrow = num_mcmc, ncol = num_repeat)
    for (irep in 1: num_repeat){
      indicator_of_coupling[ 1 : (couple_metropolis_to_condbern(alpha, sum_x, num_mcmc)$meeting_time  -1 ), irep] <- 1
    }
    tv_upper_bounds <- rowMeans(indicator_of_coupling)
    return(tv_upper_bounds)
  }
}


#' @title The Metropolis sampler for Conditional Bernoulli distribution
#' @description meeting time between Metropolis sampler in gibbs_condbern.R and the Conditional-Bernoulli distribution 
#' @param alpha vector of probabilities
#' @param sum_x Conditional-Bernoulli constraint, sum of all entries
#' @param num_repeat number of repetitions
#' @param max_iterations maximum number of iterations, at which to interrup the while loop; Inf by default
#' @return Meeting times stored in a vector
#' @export 
#' 
get_meetingtime_metropolis_idchecking <- function(alpha, sum_x, num_repeat, max_iterations = Inf){
  meeting_times <- rep(0, num_repeat)
  for (irep in 1 : num_repeat){
    meeting_time <- couple_metropolis_to_condbern(alpha = alpha, sum_x = sum_x, num_mcmc = max_iterations, compute_sse = FALSE, verbose = FALSE)$meeting_time
    meeting_times[irep] <- meeting_time 
  }
  meeting_times
}
