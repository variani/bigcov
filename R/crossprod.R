
#' @export
big_crossprod <- function(data,
  num_splits = 2,
  debug = FALSE, verbose = 0, ...)
{
  ### check
  if(class(data)[1] != "matrix") {  
    data <- as.matrix(data)
  }
  
  ### variables
  n <- nrow(data)
  p <- ncol(data)
  
  ### split into batches
  if(verbose) {
    cat(" - big_crossprod: preparing batches of block-pairs from ", num_splits, "splits\n")
  } 
  split_size <- compute_split_size(p, num_splits)
  beg <- compute_splits_beg(p, split_size)
  end <- compute_splits_end(p, split_size, beg)
  
  num_splits <- length(beg)

  batches <- llply(1:num_splits, function(i) {
    list(batch = i, beg = beg[i], end = end[i])
  })
  
  ### define the grid
  grid <-  expand.grid(1:num_splits, 1:num_splits)
  colnames(grid) <- c("cell1", "cell2")
  
  grid <- subset(grid, cell1 <= cell2)
  
  num_batches <- nrow(grid)
  
  ### computes blocks of the cov. matrix
  if(verbose) {
    cat(" - big_crossprod: computing `crossprod` by blocks:", num_batches, "batches\n")
  }
  
  blocks <- llply(1:nrow(grid), function(i) {
    if(verbose > 1) {
      cat(" - block-pair batch", i, "/", num_batches, "\n")
    }
    cell1 <- grid[i, "cell1"]
    cell2 <- grid[i, "cell2"]
    
    ind1 <- with(batches[[cell1]], seq(beg, end))
    ind2 <- with(batches[[cell2]], seq(beg, end))
    
    if(cell1 == cell2) {
      mat <- crossprod(data[, ind1])
    } else {
      mat <- crossprod(data[, ind1], data[, ind2])
    }
    
    list(cell1 = cell1, cell2 = cell2, 
      ind1 = ind1, ind2 = ind2, mat = mat)
  })
  
  ### build the output matrix block by block
 if(verbose) {
    cat(" - big_crossprod: merging the batches into", p, "x", p, "output matrix\n")
  }
  
  outmat <- matrix(0, ncol = p, nrow = p)
  
  for(i in 1:length(blocks)) {
    ind1 <- blocks[[i]]$ind1
    ind2 <- blocks[[i]]$ind2
    mat <- blocks[[i]]$mat
    
    outmat[ind1, ind2] <- mat
  }
  outmat <- Matrix::forceSymmetric(outmat, uplo = "U") 
  outmat <- as.matrix(outmat)
  
  rownames(outmat) <- colnames(data)
  colnames(outmat) <- colnames(data)
  
  ### return
  if(debug) {
   out <- list(batches = batches, grid = grid, blocks = blocks,
     mat = outmat)
  } else {
    return(outmat)
  }
}

#' @export
big_tcrossprod <- function(data,
  num_splits = 2,
  debug = FALSE, verbose = 0, ...)
{
  ### check
  if(class(data)[1] != "matrix") {  
    data <- as.matrix(data)
  }
  
  ### variables
  n <- nrow(data)
  p <- ncol(data)
  
  ### split into batches
  if(verbose) {
    cat(" - big_tcrossprod: preparing batches of block-pairs from ", num_splits, "splits\n")
  }  
  # split
  split_size <- compute_split_size(n, num_splits)
  beg <- compute_splits_beg(n, split_size)
  end <- compute_splits_end(n, split_size, beg)
  
  num_splits <- length(beg)

  batches <- llply(1:num_splits, function(i) {
    list(batch = i, beg = beg[i], end = end[i])
  })
  
  ### define the grid
  grid <-  expand.grid(1:num_splits, 1:num_splits)
  colnames(grid) <- c("cell1", "cell2")
  
  grid <- subset(grid, cell1 <= cell2)
  
  num_batches <- nrow(grid)
  
  ### computes blocks of the cov. matrix
  if(verbose) {
    cat(" - big_tcrossprod: computing `tcrossprod` by blocks:", num_batches, "batches\n")
  }
  
  blocks <- llply(1:nrow(grid), function(i) {
    if(verbose > 1) {
      cat(" - block-pair batch", i, "/", num_batches, "\n")
    }
    cell1 <- grid[i, "cell1"]
    cell2 <- grid[i, "cell2"]
    
    ind1 <- with(batches[[cell1]], seq(beg, end))
    ind2 <- with(batches[[cell2]], seq(beg, end))
    
    if(cell1 == cell2) {
      mat <- tcrossprod(data[ind1, ])
    } else {
      mat <- tcrossprod(data[ind1, ], data[ind2, ])
    }
    
    list(cell1 = cell1, cell2 = cell2, 
      ind1 = ind1, ind2 = ind2, mat = mat)
  })
  
  ### build the output matrix block by block
 if(verbose) {
    cat(" - big_tcrossprod: merging the batches into", n, "x", n, "output matrix\n")
  }
  
  outmat <- matrix(0, ncol = n, nrow = n)
  
  for(i in 1:length(blocks)) {
    ind1 <- blocks[[i]]$ind1
    ind2 <- blocks[[i]]$ind2
    mat <- blocks[[i]]$mat
    
    outmat[ind1, ind2] <- mat
  }
  outmat <- Matrix::forceSymmetric(outmat, uplo = "U") 
  outmat <- as.matrix(outmat)
  
  rownames(outmat) <- rownames(data)
  colnames(outmat) <- rownames(data)
  
  ### return
  if(debug) {
   out <- list(batches = batches, grid = grid, blocks = blocks,
     mat = outmat)
  } else {
    return(outmat)
  }
}
