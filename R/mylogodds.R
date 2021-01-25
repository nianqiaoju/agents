#' log-odds
#' @param p vector of probabilities
#' @return log(p/(1-p)) for each entry
#' @export

mylogodds <- function(p){
  log(p/(1-p))
}