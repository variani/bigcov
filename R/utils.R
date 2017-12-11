
#---------------
# Convert
#---------------

#' @export
convert_cdat <- function(cmat, diag = TRUE, lower = TRUE)
{
  stopifnot(requireNamespace("reshape2"))
  
  filter_na <- !lower
  
  ### induce NA
  if(!lower) {
    ind <- lower.tri(cmat, diag = !diag)
    cmat[ind] <- NA
  }
  
  cdat <- reshape2::melt(cmat)
  colnames(cdat) <- c("row", "col", "value")
  
  if(filter_na) {
    cdat <- subset(cdat, !is.na(value))
  }
  
  return(cdat)
}

#---------------
# Splits
#---------------

compute_split_size <- function(n, num_splits)
{
  split_size <- ceiling(n / num_splits)
  stopifnot(split_size > 1)  
  
  return(split_size)
}

compute_splits_beg <- function(n, split_size)
{
  beg <- seq(1, n, by = split_size)
  beg <- beg[beg <= n]
  
  return(beg)
}

compute_splits_end <- function(n, split_size, beg)
{
  end <- beg - 1 + split_size
  end[length(end)] <- n
  
  return(end)
}
