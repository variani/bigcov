
#' @export
bigcov <- function(data, center = TRUE, scale = FALSE,
  num_batches = 2,
  debug = FALSE, verbose = 0)
{
  ### variables
  n <- nrow(data)
  p <- ncol(data)
  
  ### scale
  data <- scale(data, center = center, scale = scale)
  
  ### split into batches
  batch_size <- ceiling(p / num_batches)
  
  beg <- seq(1, by = batch_size, length = num_batches)
  end <- beg - 1 + batch_size
  end[length(end)] <- p


  batches <- llply(1:num_batches, function(i) {
    list(batch = i, beg = beg[i], end = end[i])
  })
  
  ### define the grid
  grid <-  expand.grid(1:num_batches, 1:num_batches)
  colnames(grid) <- c("cell1", "cell2")
  
  grid <- subset(grid, cell1 <= cell2)
  
  ### computes blocks of the cov. matrix
  blocks <- llply(1:nrow(grid), function(i) {
    if(verbose) {
      cat(" - blocks pair", i, "/", nrow(grid), "\n")
    }
    cell1 <- grid[i, "cell1"]
    cell2 <- grid[i, "cell2"]
    
    ind1 <- with(batches[[cell1]], seq(beg, end))
    ind2 <- with(batches[[cell2]], seq(beg, end))
    
    if(cell1 == cell2) {
      cov <- crossprod(data[, ind1]) / (n - 1)
    } else {
      cov <- crossprod(data[, ind1], data[, ind2]) / (n - 1)
    }
    
    list(cell1 = cell1, cell2 = cell2, 
      ind1 = ind1, ind2 = ind2, cov = cov)
  })
  
  ### build the cov. matrix block by block
  covmat <- matrix(0, ncol = p, nrow = p)
  
  for(i in 1:length(blocks)) {
    ind1 <- blocks[[i]]$ind1
    ind2 <- blocks[[i]]$ind2
    cov <- blocks[[i]]$cov
    
    covmat[ind1, ind2] <- cov
    covmat[ind2, ind1] <- cov
  }
  
  rownames(covmat) <- colnames(data)
  colnames(covmat) <- colnames(data)
  
  ### return
  if(debug) {
   out <- list(batches = batches, grid = grid, blocks = blocks,
     cov = covmat)
  } else {
    return(covmat)
  }
}
