#' @title G(n,p) Erdos-Renyi model
#' @description sample a random graph with n vertices and m edges
#' @param N integer constant, the population size.
#' @param p the probability for drawing an edge between any two vertices, default to 0.5.
#' @return a matrix whose n-th row is \math{N(n)}
#' @export


## this function uses the sample_smallworld function from package igraph, but I do not completely understand the algorithms there
## according to wikipedia, I can implement my own sample_smallworld and it requires fewer parameters
network_sample_gnp <- function(N, p = 0.1){
	gnp <- igraph::sample_gnp(n = N, p = p);
	adj_gnp <- igraph::as_adjacency_matrix(gnp);
	degrees <- apply(adj_gnp, 1, sum);
	mdegree <- max(degrees);
	neighbors <- matrix(0, nrow = N, ncol = mdegree);
	for(n in 1 : N){
	  if(degrees[n] > 0)  neighbors[n, 1 : degrees[n]] <- which(adj_gnp[n, ] == 1);
	}
	neighbors
}
