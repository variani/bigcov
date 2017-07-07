
#' @export
bigdat_grm <- function(data, num_batches = NULL, batch_size = NULL,
  verbose = 0,
  ids, ...)
{
  # GRM as given in Patterson 2006, formulas (1)-(3)
  
  ### convert input `data` into an object of class `bigdat`
  bdat <- bigdat(data, num_batches = num_batches, batch_size = batch_size, ...)
  
  N <- bigdat_nrow1(bdat)
  M <- bigdat_ncols_sum(bdat)
  num_batches <- bigdat_nbatch(bdat)
  
  ### computes blocks of the cov. matrix
  if(verbose) {
    cat(" - bigdat_tcrossprod: computing `tcrossprod`:", num_batches, "batches\n")
  }
  
  product <- matrix(0, ncol = N, nrow = N)
  
  l_ply(seq(1, num_batches), function(i) {
    if(verbose > 1) {
      cat(" - batch", i, "/", num_batches, "\n")
    }
    
    mat <- bigdat_batch(bdat, batch = i)

    col_means <- colMeans(mat, na.rm = TRUE)
    col_freq <- col_means / 2
    col_sd <- sqrt(col_freq * (1 - col_freq))
    
    # center
    mat <- sweep(mat, 2, col_means, "-")
    # impute missing
    ind <- is.na(mat)
    mat[ind] <- 0
    # scale
    mat <- sweep(mat, 2, col_sd , "/")
    
    # https://stackoverflow.com/questions/2628621/how-do-you-use-scoping-assignment-in-r
    product <<- product + tcrossprod(mat)
  })
  
  product <- product / M
  
  ### ids
  if(missing(ids)) {
   ids <- as.character(1:N)
  } 
    
  rownames(product) <- ids
  colnames(product) <- ids 
  
  return(product)
}

