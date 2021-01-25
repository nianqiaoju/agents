#' Access the last element in a vector
#'
#' @param x vector 
#' @return last element in the vector x 
#'
#' @export
#' 
mylast <- function(x) x[length(x)]


#' @param x numeric vector
#' @return last no-na element in the vector, assuming that all na's are at the tail.
#' @export
mylast_nona <- function(x) x[sum(!is.na(x))]