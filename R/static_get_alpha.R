#' @title Parametrization of the static model
#' @param features N by num_features matrix
#' @param beta length is number of features
#' @return the alpha vector
#' @export  

static_get_alpha <- function(beta, features){
  if(dim(features)[1] != length(beta)) stop("features and beta not compatible")
  as.vector(1 / (1 + exp( - beta %*% features)))
}
