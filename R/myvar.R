#' @export
myvar <- function(x){
  var(x[!is.infinite(x)])
}

mysd <- function(x){
  sd(x[!is.infinite(x)])
}