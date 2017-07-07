
#' @export
bigdat_tcrossprod <- function(data, num_batches = NULL, batch_size = NULL,
  debug = FALSE, verbose = 0, ...)
{
  ### convert input `data` into an object of class `bigdat`
  bdat <- bigdat(data, num_batches = num_batches, batch_size = batch_size, ...)
  
  N <- bigdat_nrow1(bdat)
  num_batches <- bigdat_nbatch(bdat)
  
  ### computes blocks of the cov. matrix
  if(verbose) {
    cat(" - bigdat_tcrossprod: computing `tcrossprod` by blocks:", num_batches, "batches\n")
  }
  
  product <- matrix(0, ncol = N, nrow = N)
  
  l_ply(seq(1, num_batches), function(i) {
    if(verbose > 1) {
      cat(" - batch", i, "/", num_batches, "\n")
    }
    
    mat <- bigdat_batch(bdat, batch = i)
    
    # https://stackoverflow.com/questions/2628621/how-do-you-use-scoping-assignment-in-r
    product <<- product + tcrossprod(mat)
  })
  
  rownames(product) <- rownames(data)
  colnames(product) <- rownames(data)
  
  return(product)
}
