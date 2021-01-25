#' @title G(n,p) Erdos-Renyi model
#' @description sample a random graph with n vertices and m edges
#' @param N integer constant, the population size.
#' @param p the probability for drawing an edge between any two vertices, default to 0.5.
#' @export


## this function uses the sample_smallworld function from package igraph, but I do not completely understand the algorithms there
## according to wikipedia, I can implement my own sample_smallworld and it requires fewer parameters
network_sample_gnp <- function(N, p = 0.5){
	gnp <- igraph::sample_gnp(n = N, p = p);
	adj_gnp <- igraph::as_adjacency_matrix(gnp);
	diag(adj_gnp) <- rep(1,N);
	adj_gnp <- apply(adj_gnp,2,function(col_) col_ / sum(col_));
	return(adj_gnp);
}
