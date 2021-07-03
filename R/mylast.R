#' @title Access the last element in a vector
#' @param x vector 
#' @param rm.na binary, if TRUE, then ignore the NA values in x.
#' @return last element in the vector x 
#' 
mylast <- function(x, rm.na = FALSE){
  if(rm.na){ ## if ignore NA
    ## this is assuming that all NA's are at the end of the vector.
    return(x[sum(!is.na(x))])
  }else{
    return(x[length(x)]);
  }
}
