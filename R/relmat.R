
#' @export
big_relmat <- function(data, center = TRUE, scale = FALSE,
  num_splits = 2, verbose = 0,
  ids, ...)
{
  ### prepare the matrix of genotypes: to be centered / scaled
  if(class(data)[1] != "matrix") {
    data <- as.matrix(data)
  }
  
  data <- scale(data, center = center, scale = scale)

  ### var
  M <- ncol(data)
  N <- nrow(data)
  
  ### compute the var-covar matrix
  relmat <- big_tcrossprod(data, num_splits = num_splits, verbose = verbose, ...)
  relmat <- relmat / M
  
  ### ids
  if(missing(ids)) {
   ids <- as.character(1:N)
  } 
    
  rownames(relmat) <- ids
  colnames(relmat) <- ids 
  
  return(relmat)
}
