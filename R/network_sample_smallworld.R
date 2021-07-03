#' @title Small world graph
#' @description sample a small world graph and return the `neighbors' matrix representing the graph.
#' @param N integer constant, the population size.
#' @param dim integer constant, dimension of the starting lattice;
#' @param nei integer constant, the neighborhood within which the vertices of the lattice will be connected.
#' @param p real constant between zero and one, the rewiring probability.
#' @return a matrix whose n-th row is \math{N(n)}
#' @export


## this function uses the sample_smallworld function from package igraph, but I do not completely understand the algorithms there
## according to wikipedia, I can implement my own sample_smallworld and it requires fewer parameters
network_sample_smallworld <- function(N, dim = 1, nei = 1, p = 0.3){
	sw <- igraph::sample_smallworld(dim = dim, size = N, nei = nei, p = p);
	adj_sw <- igraph::as_adjacency_matrix(sw);
	degrees <- apply(adj_sw, 1, sum);
	mdegree <- max(degrees);
	neighbors <- matrix(0, nrow = N, ncol = mdegree);
	for(n in 1 : N){
	  if(degrees[n] > 0) neighbors[n, 1 : degrees[n]] <- which(adj_sw[n, ] == 1);
	}
	neighbors
}
