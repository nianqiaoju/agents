#' @title Twisted SIR kernels
#' @description  samples from xt given it and x(t-1). This function is used in sir_apf.R.
#' @param xxprev agent states at the previous time;
#' @param it infection counts at the current time;
#' @param model_config a list containing model parameters, must pass the test of check_model_config;
#' @return a sample of xxprev, vector of length N.
#' @export

sir_kernel_twisted_i <- function(xxprev, it, model_config){
  xx <- rep(0,length(xxprev))
  xx[which(xxprev == 2)] <- 2 # the recovered stay recovered
  xx[which(xxprev != 2)] <- rcondbern(sum_x = it, alpha = sir_get_alpha_i(xxprev, model_config)[which(xxprev != 2)])
  xx[which(xxprev == 1)] <- 2 - xx[which(xxprev ==1)] # if an infected turns not infected, it must have recovered.
  return(xx)
}