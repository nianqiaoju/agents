#' @title Small world graph
#' @description sample a small world graph and return the (normalized) adjacency matrix
#' @param N integer constant, the population size.
#' @param dim integer constant, dimension of the starting lattice;
#' @param nei integer constant, the neighborhood within which the vertices of the lattice will be connected.
#' @param p real constant between zero and one, the rewiring probability.
#' @export


## this function uses the sample_smallworld function from package igraph, but I do not completely understand the algorithms there
## according to wikipedia, I can implement my own sample_smallworld and it requires fewer parameters
network_sample_smallworld <- function(N, dim = 1, nei = 1, p = 0.3){
	sw <- igraph::sample_smallworld(dim = dim, size = N, nei = nei, p = p);
	adj_sw <- igraph::as_adjacency_matrix(sw);
	diag(adj_sw) <- rep(1,N);
	adj_sw <- apply(adj_sw,2,function(col_) col_ / sum(col_));
	return(adj_sw);
}
