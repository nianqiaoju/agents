## 
static_xx_gibbs <- function(xx, alpha, rho, y){
  xxnotn <- sum(xx[-1])
  for (n in 1 : N){
    # check if xx[n] has to be 1
    if (xxnotn < y){
      xx[n] <- 1
    }else{
      p1 <- alpha[n] * dbinom(x = y, size = 1 + xxnotn, prob = rho)
      p0 <- (1 - alpha[n]) * dbinom(x = y, size = xxnotn, prob = rho)
      xx[n] <- (runif(1) < ( p1 / (p1 + p0)))
      # if (n == 1) print( p1 / (p1 + p0)) ## debug print message
    }
    xxnotn <- xxnotn + xx[n] - xx[min(N, n+1)]
  }
  return(xx)
  # static_xx_gibbs_cpp(xx, alpha, rho, y)
}


static_xx_blocked_gibbs <- function(xx, alpha, rho, y, block_size = 20){
  blocks <- split(1 : length(xx), ceiling(seq_along(1 : length(xx))/block_size))
  num_blocks <- length(blocks)
  for(iblock in 1 : num_blocks){
    sum_rest <- sum(xx[-blocks[[iblock]]])
    current_support <- 0 : length(blocks[[iblock]])
    logv <- logdpoisbinom_cpp(alpha = alpha[blocks[[iblock]]]) + dbinom(x = y, size = current_support + sum_rest, prob =  rho, log = TRUE)
    sum_inblock <- sample(current_support, size = 1, replace = FALSE, prob = lw.normalize(logv))
    xx[blocks[[iblock]]] <- rcondbern(sum_x = sum_inblock, alpha = alpha[blocks[[iblock]]])
  }
  return(xx)
}
