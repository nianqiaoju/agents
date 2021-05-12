#' @title Coupled Metropolis algorithm for the Conditional-Bernoulli distribution
#' @description A coupling of two metropolis samplers targetting the same Conditional-Bernoulli distribution
#' Chain1 starts from the stationary distribution and chain2 is initialized uniformly from the support
#' @param alpha vector of probabilities for the Poisson-Binomial distribution 
#' @param sum_x sum of x
#' @param num_mcmc number of MCMC iterations; if num.mcmc = Inf, then run the algorithm until the two chains couple
#' @param compute_sse binary, default is FALSE; if compute_sse = TRUE, then returns the sum-squared-error distance from chain1 to chain2
#' @param verbose binary, whether to print messages
#' @return a list that contains 
#' \itemize{
#' \item meeting time 
#' \item sum-squared-error distance from chain1 to chain2
#' }
#' @export

couple_metropolis_to_condbern <- function(alpha, sum_x, num_mcmc, compute_sse = FALSE, verbose = FALSE){
  N <- length(alpha)
  coupled <- FALSE
  meeting_time <- Inf
  compute_see <- (is.finite(num_mcmc) & compute_sse)
  ## if num.mcmc is not Inf, then save the vector of sum-squared-erros
  if(compute_sse){
    sse_vector <- rep(0, num_mcmc)
  }else{
    sse_vector <- NA 
  }
  ## simulate chain 1 from stationarity
  xx1 <- idchecking_cpp(sum_x,alpha, runif(N));
  ## simulate chain 2 from random initialization
  xx2 <- rep(0, N)
  xx2[sample.int(n = N, size = sum_x, replace = FALSE)] <- 1
  iter <- 1 
  while(iter <= num_mcmc & !coupled){
    one_iteration <- couple_condbern_kernel(xx1,xx2,sum_x,sum_x,alpha,alpha)
    xx1 <- one_iteration[1,]
    xx2 <- one_iteration[2,]
    if(compute_sse){
      sse_vector[iter] <- sum((xx1 - xx2)**2)
    }
    if (all(xx1 == xx2)){ ### it couples here
      if (verbose)  cat("coupled at iteration ",iter,'\n')
      coupled <- TRUE
      meeting_time <- iter
      break
    }
    iter <- iter + 1
  }
  return(list(sse_vector = sse_vector, meeting_time = meeting_time))
}
