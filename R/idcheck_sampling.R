#' @title ID-check sampling
#' @description  obtain a sample from the conditional Bernoulli distribution 
#' uses the recursive method 2 from Chen and Liu (1997)
#' @param sum_x sum of the samples 
#' @param alpha vector of probabilities 
#' @param config a list containing
#' \itemize{
#' \item num_sample : number of samples from conditional-bernoulli; default is 1
#' }
#' @return a matrix or a vector containing the Conditional Bernoulli samples
#' @export
#' 

idcheck_sampling <- function(sum_x, alpha, config = list(num_sample = 1)){
  N_ <- length(alpha)
  if(sum_x == N_) return(matrix(1, nrow = config$num_sample, ncol = N_))
  if(sum_x == 0) return(matrix(0, nrow = config$num_sample, ncol = N_))
  ## compute all the intermediate q(h,n) in log scale 
  ## q(h,n) = PB(h; alpha[n : N])
  ## base cases: 
  ## (1) q(0,n) = prod(1 - alpha[n:N]) for all n 
  ## (2) q(h, N - h + 1) = prod(alpha[N - h + 1 : N])
  ## (3) q(h, n ) = 0 for h > N - n + 1
  ## recursion: q(h,n) = alpha[n] * q(h-1, n+1) + (1 - alpha[n]) * q(h, n+1)
  ## logq[h + 1, n] = log(q(h, n)) because the row index is off-set by 1
  logq <- matrix(-Inf, nrow = sum_x + 1 , ncol = N_)
  logq[1 , ] <- rev(cumsum(log(1 - rev(alpha)))) # base case (1)
  logq[1 + 1 , N_] <- log(alpha[N_]) # base case (2) for h = 1
  for (n in (N_ - 1) : 1){
    for (h in 1 : min(sum_x, N_ - n + 1)){
      ls1 <- log(alpha[n]) + logq[h - 1 + 1, n + 1]
      ls2 <- log(1 - alpha[n]) + logq[h + 1, n + 1]
      maxls <- max(ls1, ls2)
      if(is.infinite(maxls)){
        logq[h + 1, n] <- - Inf 
      }else{
        logq[h + 1, n] <- log(exp(ls1 - maxls) + exp(ls2 - maxls)) + maxls
      }
    }
  }
  ## can check if logq[sum_x + 1, 1] = dpoisbinom(sum_x, alpha) or any other density [in test/test_idcheck_sampling.R]
  ## start to sample from conditional Bernoulli
  xx <- rep(0,N_) ## initialize the sample
  h <- 0 ## partial sum 
  ## probablity of xx[n] = 1 is alpha[n] * q(i - h - 1, n + 1) / q(i - h, n)
  for (n in 1: (N_ - 1)){
    logp <- log(alpha[n]) + logq[sum_x - h - 1 + 1, n + 1] - logq[sum_x - h + 1, n]
    if (log(runif(1)) < logp){  
      xx[n] <- 1
      h <- h + 1
    }
    if (h == sum_x){# exit for-loop whenever the partial sum is the total sum
      break
    }
  }
  ## final step is deterministic
  ## we should have partial_sum = sum_x - 1 or sum_x and xx[N_] = 1 or 0;
  ## otherwise there is a problem
  if (h < sum_x - 1 | h > sum_x){
    warning("total sum cannot be sum_x")
  }
  if (h == sum_x - 1){
    xx[N_] <-  1
  }
  if (config$num_sample == 1){
    return(xx)
  }else{
    xx_matrix <- matrix(0, nrow = config$num_sample, ncol = N_)
    xx_matrix[1, ] <- xx 
    for (isample in 2 : config$num_sample){
      h <- 0 ## partial sum 
      ## probablity of xx[n] = 1 is alpha[n] * q(i - h - 1, n + 1) / q(i - h, n)
      for (n in 1: (N_ - 1)){
        # cat("index = ", ii, "current partial sum = ", r, "\n")
        logp <- log(alpha[n]) + logq[sum_x - h - 1 + 1, n + 1] - logq[sum_x - h + 1, n]
        if (log(runif(1)) < logp){  
          xx_matrix[isample,n] <- 1
          h <- h + 1
        }
        if (h == sum_x){
          break
        }
      }
      if (h < sum_x - 1 | h > sum_x){
        warning("total sum cannot be sum_x")
      }
      if (h == sum_x - 1){
        xx_matrix[isample, N_] <- 1
      }
    }
    return(xx_matrix)
  }
}
