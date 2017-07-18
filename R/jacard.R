
#' @param binary Data are binary attributes or not? Default, \code{FALSE}.
#'    The binary format of data is assumed to compute the Jacard similarity matrix.
#'    In general case, all non-zero values are converted to binaries entries.
#'    If your data is binary, set to \code{TRUE}. 
#'
#' @export
jacard <- function(data, ids, sparse = TRUE, binary = FALSE, ...)
{
  # Cortesy of https://stats.stackexchange.com/a/89947
  # See also https://stats.stackexchange.com/a/61910
    
  # convert data to objects of Matrix class, 
  # as computation is a way faster for sparse data
  M <- Matrix(data, sparse = sparse)
  
  # convert to binary from any non-zero values
  if(!binary) {
    M@x <- rep(1, length(M@x))
  }
  
  ### main computation based on code at https://stats.stackexchange.com/a/89947
  A <- tcrossprod(M)
  
  # pairs with non-zero common values
  ind <- which(A > 0, arr.ind = TRUE)
  ind1 <- ind[, 1]
  ind2 <- ind[, 2]
  
  sums <- rowSums(M)
  
  # common values
  vals <- A[ind]
  
  # Jacard formula: #common / (#i + #j - #common)
  J <- sparseMatrix(i = ind1, j = ind2,
    x = vals / (sums[ind1] + sums[ind2] - vals),
    dims = dim(A))
  
  as.matrix(J)
}

#' @export
bigdat_jacard <- function(data, num_batches = NULL, batch_size = NULL,
  check_na = FALSE,
  sparse = TRUE,
  verbose = 0,
  ids, ...)
{
  stop("depreciated")
  
  ### args
  stopifnot(!check_na)
  
  ### convert input `data` into an object of class `bigdat`
  bdat <- bigdat(data, num_batches = num_batches, batch_size = batch_size, ...)
  
  N <- bigdat_nrow1(bdat)
  M <- bigdat_ncols_sum(bdat)
  num_batches <- bigdat_nbatch(bdat)
  
  ### map: compute blocks of the cov. matrix
  if(verbose) {
    cat(" - bigdat_jacard: computing:", num_batches, "batches\n")
  }
  
  out <- llply(seq(1, num_batches), function(i) {
    if(verbose > 1) {
      cat(" - batch", i, "/", num_batches, "\n")
    }
    
    # get a batch of data
    mat <- bigdat_batch(bdat, batch = i)

    M <- Matrix(mat, sparse = sparse)

    list(product = tcrossprod(M), sums = rowSums(M))
  })
  
  ### reduce
  if(length(out) == 1) {
    product <- out[[1]][["product"]]
    sums <- out[[1]][["sums"]]
  } else {
    product <- out[[1]][["product"]]
    sums <- out[[1]][["sums"]]
    for(i in seq(2, length(out))) {
      product <- product + out[[i]][["product"]]
      sums <- sums + out[[i]][["sums"]]
    }
  }
  
  ### compute Jacard
  # pairs with non-zero common values
  ind <- which(product > 0, arr.ind = TRUE)
  ind1 <- ind[, 1]
  ind2 <- ind[, 2]
  
  # common values
  vals <- product[ind]
  
  # Jacard formula: #common / (#i + #j - #common)
  J <- sparseMatrix(i = ind1, j = ind2,
    x = vals / (sums[ind1] + sums[ind2] - vals),
    dims = dim(product))
  
  ### covert back to matrix  
  product <- as.matrix(J)
  
  ### ids
  if(missing(ids)) {
   ids <- as.character(1:N)
  } 
    
  rownames(product) <- ids
  colnames(product) <- ids 
  
  return(product)
}

