#' @export
bigcov <- function(data, center = TRUE, scale = FALSE,
  num_splits = 2,
  verbose = 0, ...)
{
  big_cov(data, center = center, scale = scale,
  num_splits = num_splits,
  verbose = verbose, ...)
}

#' @export
big_cov <- function(data, center = TRUE, scale = FALSE,
  num_splits = 2,
  verbose = 0, ...)
{
  ### variables
  n <- nrow(data)
  p <- ncol(data)
  
  transform <- center | scale
    
  ### scale
  if(transform) {
    data <- scale(data, center = center, scale = scale)
  }
  
  ### compute crossprod
  covmat <- big_crossprod(data,
    num_splits = num_splits, verbose = verbose, ...)
  
  covmat <- covmat / (n - 1)
 
  ### return
  return(covmat)
}
