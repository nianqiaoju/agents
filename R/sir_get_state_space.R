#' @title State space for SIR model
#' @description create the space \{0,1,2\}^N
#' @param N population size
#' @return the space, as a 3**N by N matrix 
#' @examples sir_get_state_space(3)
#' @export

sir_get_state_space <- function(N){
  # returns the space {0,1}^N
  state_space <- as.matrix(expand.grid(rep(list(0L:2L), N)))
  colnames(state_space) <- c(1:N)
  return(state_space)
}
