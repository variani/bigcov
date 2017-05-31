
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
