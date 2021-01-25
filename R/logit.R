#' @title Logit function 
#' @export
logit <- function(x){
  return(log(x/(1-x)))
}

#' @title Expit
#' @description  expit function, inverse of logit 
#' 
#' @export
expit <- function(z) 1 / (1 + exp(-z));
