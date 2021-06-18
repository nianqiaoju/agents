#' @title Normalize log weights using the log-sum-exp trick.
#' @param lw log weights
#' @return \code{lw.normalize} returns normalized weights.
#' @export

lw.normalize <- function(lw){
  # output: normalized weights, aka they must sum to 1 
  # if all the weights are zero, then raise an error message.
  maxlw <- max(lw);
  if(is.infinite(maxlw)){
    warning("all the weights are zero!");
  }else{
    lw <- lw - maxlw;
    return(exp(lw)/sum(exp(lw)));
  }
}

#' @rdname lw.normalize
#' @return \code{lw.normalize} computes the log normalizing constant 
#' @export
lw.logsum <- function(lw){
  # input: lw - log weights 
  # output: sum of weights, i.e. log(sum(exp(lw)))
  maxlw <- max(lw)
  if(is.infinite(maxlw)){
    # warning("all the weights are zero");
    return(-Inf);
  }
  else{
    return(maxlw + log(sum(exp(lw - maxlw))));
  }
}

#' @rdname lw.normalize
#' @return \code{lw.logmean} computes log of mean of the vector \eqn{exp(lw)}.
#' @export
#' 
lw.logmean <- function(lw){
  lw.logsum(lw) - log(length(lw));
}

#' @rdname lw.normalize
#' @return \code{lw.logsd} computes log of standard deviation of the vector \eqn{exp(lw)}.
#' @export
lw.logsd <- function(lw){
  0.5 * lw.logmean(2 * log(abs(exp(lw) - exp(lw.logmean(lw)))));
}