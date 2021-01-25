#' @title Conditional Bernoulli distribution 
#' @description Gives one sample from the conditional Bernoulli distribution, either with exact method or with MCMC. 
#' This function is a wrapper for idchecking_cpp and metropolis_condbern_cpp
#' @param sum_x sum of the Bernoullis
#' @param alpha vector of success probabilities
#' @param exact binary, if TRUE, then use id-check sampling, otherwise use the metropolis algorithm. Default is TRUE.
#' @param num_sample default is NULL. 
#' @return sample, a length N binary vector, with 1 at the chosen entry
#'
#' @export
#' 
rcondbern <- function(sum_x, alpha, exact = TRUE, num_sample = NULL){
  if (sum_x == length(alpha)) return(rep(T, sum_x))
  if (sum_x == 0) return(rep(F, length(alpha)))
  if (exact){
    return(idchecking_cpp(sum_x = sum_x, alpha = alpha, random_number = runif(length(alpha))));
  }else{
    if (is.null(num_sample)){
      num_sample <- ceiling(length(alpha) * log(length(alpha))); ## can tune this parameter
    }
    return(metropolis_condbern_cpp(sum_x, alpha, num_sample))
  }
}
