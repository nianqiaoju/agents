#' Normalize log weights
#'
#' @param lw log weights
#' @return normalized weights
#' @export

lw.normalize <- function(lw){
  # input: lw - log weights 
  # output: normalized weights, aka they must sum to 1 
  maxlw <- max(lw)
  if(is.infinite(maxlw)){
    return(rep(1 / length(lw), length(lw))); 
    ## if all the weights are 0, returnrs rep(1/n,n). this might be very dangerous!
  }else{
    lw <- lw - maxlw
    return(exp(lw)/sum(exp(lw)))
  }
  # return(exp(lw)/exp(lw.logsum(lw)))
}

#' log-normalizing constant
#' @param lw log-weights
#' @return log-normalizing constant
#' @export
lw.logsum <- function(lw){
  # input: lw - log weights 
  # output: sum of weights, i.e. log(sum(exp(lw)))
  maxlw <- max(lw)
  if(is.infinite(maxlw)){
    return(-Inf)
  }
  else{
    return(maxlw + log(sum(exp(lw - maxlw))))
  }
}

#' log-mean 
#' @param lw log-weights
#' @return log of the mean 
#' @export
#' 
lw.logmean <- function(lw){
  lw.logsum(lw) - log(length(lw))
}

#' log-sd
#' @param  lw log-weights
#' @return log of the standard deviation 
#' @export
lw.logsd <- function(lw){
  0.5 * lw.logmean(2 * log(abs(exp(lw) - exp(lw.logmean(lw)) )))
}