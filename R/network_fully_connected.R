#' @title Network structure
#' @description generate the (normalized) adjacency matrix of a fully connected network
#' @param N integer constant, the population size.
#' @export
network_fully_connected <- function(N){
  adjacency_matrix_a <- matrix(1, nrow = N, ncol = N)
  # diag(adjacency_matrix_a) <- 0 ## if a agent cannot be its own neighbor
  adjacency_matrix_b <- apply(adjacency_matrix_a, 2, function(col_) col_ / sum(col_))
  adjacency_matrix_b
}


