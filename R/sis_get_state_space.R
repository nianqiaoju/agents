#' State space of SIS model
#' @description  Please only use this function for N <= 5.
#' @param N population size 
#' @return the state spacespace \{0,1\}^N
#' @examples all_x_sis(3)
#' @export


sis_get_state_space <- function(N){
  # returns the space {0,1}^N
  M <- as.matrix(expand.grid(rep(list(0L:1L), N)))
  colnames(M) <- c(1:N)
  return(apply(M, c(1,2), as.logical))
}
