
#' @export
bigdat_grm <- function(data, grm = c("Patterson", "Jacard"),
  num_batches = NULL, batch_size = NULL,
  filter = TRUE, 
  check_na = TRUE, min_callrate = 0.98, filter_mono = TRUE, maf_min = NULL, maf_max = NULL, 
  sparse = TRUE, 
  cores = 1, 
  # debugging
  num_batches_exec = NULL,
  verbose = 0,
  ids, ...)
{
  # GRM as given in Patterson 2006, formulas (1)-(3)
  
  ### args
  grm <- match.arg(grm)
  
  ### vars
  parallel <- (cores > 1)

  if(filter) {
    filter_callrate <- !is.null(min_callrate)
    filter_maf_min <- !is.null(maf_min)
    filter_maf_max <- !is.null(maf_max)
    filter_maf <- filter_maf_min | filter_maf_max

    if(verbose > 0) {
      cat("  -- filters:", ifelse(filter_callrate, "callrate", ""), 
        "/", ifelse(filter_mono, "mono", ""), 
        "/", ifelse(filter_maf_min, "maf_min", ""), 
        "/", ifelse(filter_maf_max, "maf_max", ""), 
        "/", ifelse(check_na, "check_na", ""), 
        "\n")  
    }
  }
  
  ### convert input `data` into an object of class `bigdat`
  bdat <- bigdat(data, num_batches = num_batches, batch_size = batch_size, ...)
  
  N <- bigdat_nrow1(bdat)
  M <- bigdat_ncols_sum(bdat)
  num_batches <- bigdat_nbatch(bdat)
  
  ### map: compute blocks of the cov. matrix
  if(verbose) {
    cat(" - bigdat_tcrossprod: computing `tcrossprod`:", num_batches, "batches\n")
  }
  
  # output matrix shared in memory 
  product <- big.matrix(nrow = N, ncol = N, init = 0.0, type = 'double')
  mutex_product <- boost.mutex()

  if(parallel) {
    registerDoParallel(cores = cores)
  }
  
  num_batches_loop <- ifelse(is.null(num_batches_exec), num_batches, num_batches_exec)
  stopifnot(num_batches_loop <= num_batches)
  
  out <- llply(seq(1, num_batches_loop), function(i) {
    if(verbose > 1) {
      cat(" - batch", i, "/", num_batches, "\n")
    }
    
    # get a batch of data
    mat <- bigdat_batch(bdat, batch = i)
    
    # estimate freqs
    col_means <- colMeans(mat, na.rm = TRUE)
    col_freq <- col_means / 2
    col_sd <- sqrt(col_freq * (1 - col_freq))
    
    # prepare num_*
    num_markers <- ncol(mat)
    
    # filter
    num_markers_clean <- num_markers
    
    if(filter) {
      num_filtered_callrate <- 0
      num_filtered_all_na <- 0
      num_filtered_mono <- 0
      num_filtered_maf <- 0
      num_filtered_na <- 0
    
      if(num_markers_clean > 0 & filter_callrate) {
        col_callrate <- apply(mat, 2, function(x) sum(!is.na(x))) / nrow(mat)
      
        ind_out <- (col_callrate < min_callrate)

        if(any(ind_out)) {
          num_filtered_callrate <- sum(ind_out)
          num_markers_clean <- num_markers_clean - num_filtered_callrate
        
          if(verbose > 1) {
            cat("   --- filtered callrate markers:", num_filtered_callrate, "/", num_markers, "\n")
          }
        
          ind_in <- !ind_out
  
          mat <- mat[, ind_in]
          col_means <- col_means[ind_in]
          col_freq <- col_freq[ind_in]
          col_sd <- col_sd[ind_in]
        }      
      }
    
      if(num_markers_clean > 0 & !filter_callrate) {
        ind_out <-  is.na(col_means)
      
        if(any(ind_out)) {
          num_filtered_all_na <- sum(ind_out)
          num_markers_clean <- num_markers_clean - num_filtered_all_na
        
          if(verbose > 1) {
            cat("   --- filtered all na markers:", num_filtered_all_na, "/", num_markers, "\n")
          }
        
          ind_in <- !ind_out

          mat <- mat[, ind_in]
          col_means <- col_means[ind_in]
          col_freq <- col_freq[ind_in]
          col_sd <- col_sd[ind_in]
        }
      }
    
     if(num_markers_clean > 0 & filter_mono & !filter_maf_min) {
        ind_out <- (col_means == 0)
      
        if(any(ind_out)) {
          num_filtered_mono <- sum(ind_out)
          num_markers_clean <- num_markers_clean - num_filtered_mono
        
          if(verbose > 1) {
            cat("   --- filtered mono. markers:", num_filtered_mono, "/", num_markers, "\n")
          }
        
          ind_in <- !ind_out

          mat <- mat[, ind_in]
          col_means <- col_means[ind_in]
          col_freq <- col_freq[ind_in]
          col_sd <- col_sd[ind_in]
        }
      } 
    
      if(num_markers_clean > 0 & filter_maf) {
        if(filter_maf_min & filter_maf_max) {
          ind_out <- (col_freq < maf_min | col_freq > maf_max)
        } else if(filter_maf_min) {
          ind_out <- (col_freq < maf_min)
        } else if(filter_maf_max) {
          ind_out <- (col_freq > maf_max)
        } else {
          stop("filter_maf")
        }
      
       if(any(ind_out)) {
          num_filtered_maf <- sum(ind_out)
          num_markers_clean <- num_markers_clean - num_filtered_maf
        
          if(verbose > 1) {
            cat("   --- filtered by maf markers:", num_filtered_maf, "/", num_markers, "\n")
          }
        
          ind_in <- !ind_out
  
          mat <- mat[, ind_in]
          col_means <- col_means[ind_in]
          col_freq <- col_freq[ind_in]
          col_sd <- col_sd[ind_in]
        }
      }
    } # end of `if(filter)`
    
    ### compute GRM/Jacard
    if(verbose > 1) {
      cat("  -- clean markers ready for the analysis:", num_markers_clean, "/", num_markers, "\n") 
    }
    row_sums <- NULL  
    
    if(grm %in% c("Patterson")) {
      # center
      mat <- sweep(mat, 2, col_means, "-")
      # impute missing
      if(check_na) {
        ind <- is.na(mat)
        mat[ind] <- 0
        
        if(verbose > 1) {
          cat("  -- # NAs", sum(is.na(mat)), "\n")
        }
      }
      # scale
      mat <- sweep(mat, 2, col_sd , "/")
    
      # write result into `product`
      #product_batch <- tcrossprod(mat)
      lock(mutex_product)
      product[, ] <- product[, ] + tcrossprod(mat)
      unlock(mutex_product)
    } else if(grm %in% "Jacard") {
       if(check_na) {
        ind_out <- apply(mat, 2, function(x) any(is.na(x)))
        
        if(any(ind_out)) {
          num_filtered_na <- sum(ind_out)
          if(verbose > 1) {
            cat("   --- filtered by na markers:", num_filtered_na, "/", num_markers, "\n")
          }
        
          ind_in <- !ind_out

          mat <- mat[, ind_in]
          col_means <- col_means[ind_in]
          col_freq <- col_freq[ind_in]
          col_sd <- col_sd[ind_in]
        }
      }
      
      if(sparse) {
        mat <- Matrix(mat, sparse = TRUE)
      }
      
      # write result into `product`
      lock(mutex_product)
      product[, ] <- product[, ] + tcrossprod(mat)
      unlock(mutex_product)
            
      row_sums <- rowSums(mat)
    } else {
      stop("grm")
    }  
    
    list(batch = i, 
      num_markers = num_markers, num_markers_clean = num_markers_clean,
      row_sums = row_sums)
  }, .parallel = parallel)
  
  if(parallel) {
    ret <- stopImplicitCluster()
  }
  
  ### extract statistics from batches
  num_markers <- sapply(out, function(x) x$num_markers) %>% sum
  num_markers_clean <- sapply(out, function(x) x$num_markers_clean) %>% sum  

  if(verbose > 0) {
    cat(" - clean markers used in the analysis:", num_markers_clean, "/", num_markers, "\n") 
  }
      
  ### process output from batches
  if(grm %in% c("Patterson")) {
    product <- as.matrix(product)
    product <- product / M
  } else if(grm %in% "Jacard") {
    # `row_sums`
    row_sums <- out[[1]][["row_sums"]]
    for(i in seq(2, length(out))) {
      row_sums <- row_sums + out[[i]][["row_sums"]]
    }
    
    # covert back to matrix  
    product <- as.matrix(product)
        
    # pairs with non-zero common values
    ind <- which(product > 0, arr.ind = TRUE)
    ind1 <- ind[, 1]
    ind2 <- ind[, 2]
  
    # common values
    vals <- product[ind]
    
    # Jacard formula: #common / (#i + #j - #common)
    J <- sparseMatrix(i = ind1, j = ind2,
      x = vals / (row_sums[ind1] + row_sums[ind2] - vals),
      dims = dim(product))
  
    # covert back to matrix  
    product <- as.matrix(J)
  } else {
    stop("grm")
  }
          
  ### ids
  if(missing(ids)) {
   ids <- as.character(1:N)
  } 
    
  rownames(product) <- ids
  colnames(product) <- ids 
  
  return(product)
}

