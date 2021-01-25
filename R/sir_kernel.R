#' transition kernel for the SIR process 
#' @description This is the one step transition kernel
#' @param  xxprev which is $X_{t-1}$, length N vector presenting the current states 
#' @param  model_config a list containing model paramters
#' @return xx which is $X_t$
#' @export

sir_kernel <- function(xxprev, model_config){
  ##  need the function state_to_alpha
  alpha_t <- sir_get_alpha_plusone(xxprev, model_config);
  return(xxprev + (runif(N) < alpha_t));
}
