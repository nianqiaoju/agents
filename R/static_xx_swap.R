#' @title Swap move for Conditional Bernoulli target
#' @description For a vector of length N, performs N swap moves. This algorithm runs in O(N) operations.
#' @param xx current sample
#' @param alpha vector of success probabilities
#' @return xx_new
#' @export 

static_xx_swap <- function(xx, alpha){
  ## in these cases the problem is not interseting
	if(all(xx == 1) | all(xx == 0)) return(xx);
  ## find N and I
  N_ <- length(xx);
  sum_xx <- sum(xx);
  ## convert alpha to odds
  w <- log(alpha / (1 - alpha));
  ## convert xx to sets
  s0 <- which(xx == 0); ## these are the locations for which xx[i] = 0 and this is a vector of length N - I;
  s1 <- which(xx == 1); ## these are the locations for which xx[i] = 1 and this is a vector of length I;
  for (iter in 1 : N_){ ## N iterations
    ## choose i0 and i1
    i0_loc <- 1 + floor(runif(1) * (N_ - sum_xx));
    i1_loc <- 1 + floor(runif(1) * (sum_xx));
    ## compute probability
    accept_probabilty <- w[s0[i0_loc]] / w[s1[i1_loc]] ; 
    if (runif(1) < accept_probabilty){
      ## switch indices
      i0 <- s0[i0_loc]; ## temp storage, since R syntax does not support simultaneous swaps
      s0[i0_loc] <- s1[i1_loc];
      s1[i1_loc] <- i0;
    }
  }
  ## convert s0 and s1 to binary vector representation
  xx[s0] <- 0;
  xx[s1] <- 1;
	return(xx)
}