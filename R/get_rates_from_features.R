#' @title Parameters for the SIS and SIR process
#' @param features N by num_features matrix
#' @param beta a vector of coefficients. Its length equals the number of features.
#' @return a vector expit(beta * features)
#' @export  

get_rates_from_features <- function(beta, features){
  if(dim(features)[1] != length(beta)) stop("features and beta not compatible")
  as.vector(1 / (1 + exp( - beta %*% features)));
}
