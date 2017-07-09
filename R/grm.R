
#' @export
bigdat_grm <- function(data, num_batches = NULL, batch_size = NULL,
  check_na = TRUE, filter_mono = TRUE, maf_min = NULL, 
  verbose = 0,
  ids, ...)
{
  # GRM as given in Patterson 2006, formulas (1)-(3)
  
  ### vars
  filter_maf <- !is.null(maf_min)
  stopifnot(!filter_maf)
  
  ### convert input `data` into an object of class `bigdat`
  bdat <- bigdat(data, num_batches = num_batches, batch_size = batch_size, ...)
  
  N <- bigdat_nrow1(bdat)
  M <- bigdat_ncols_sum(bdat)
  num_batches <- bigdat_nbatch(bdat)
  
  ### map: compute blocks of the cov. matrix
  if(verbose) {
    cat(" - bigdat_tcrossprod: computing `tcrossprod`:", num_batches, "batches\n")
  }
  
  out <- llply(seq(1, num_batches), function(i) {
    if(verbose > 1) {
      cat(" - batch", i, "/", num_batches, "\n")
    }
    
    # get a batch of data
    mat <- bigdat_batch(bdat, batch = i)
    
    # estimate freqs
    col_means <- colMeans(mat, na.rm = TRUE)
    col_freq <- col_means / 2
    col_sd <- sqrt(col_freq * (1 - col_freq))
    
    # filter
    if(filter_mono & !filter_maf) {
      ind_out <- (col_means == 0)
      if(any(ind_out)) {
        if(verbose > 1) {
          cat("  -- filtered mono. markers:", sum(ind_out), "\n")
        }
        
        ind_in <- !ind_out

        mat <- mat[, ind_in]
        col_means <- col_means[ind_in]
        col_freq <- col_freq[ind_in]
        col_sd <- col_sd[ind_in]
      }
    }
    
    # center
    mat <- sweep(mat, 2, col_means, "-")
    # impute missing
    if(check_na) {
      ind <- is.na(mat)
      mat[ind] <- 0
    }
    # scale
    mat <- sweep(mat, 2, col_sd , "/")
    
    if(verbose > 1) {
      cat("  -- # NAs", sum(is.na(mat)), "\n")
    }
    
    tcrossprod(mat)
  })
  
  ### reduce
  if(length(out) == 1) {
    product <- out[[1]]
  } else {
    product <- Reduce("+", out)
  }
  product <- product / M
  
  ### ids
  if(missing(ids)) {
   ids <- as.character(1:N)
  } 
    
  rownames(product) <- ids
  colnames(product) <- ids 
  
  return(product)
}

