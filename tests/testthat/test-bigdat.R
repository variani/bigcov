context("bigdat")

test_that("basic example (matrix)", {
  nrows <- 10
  ncols <- 5
  
  mat <- list(matrix(1, nrows, ncols), matrix(2, nrows, ncols))
  bdat <- bigdat(mat, batch_size = 2)
  
  for(i in seq(1, bigdat_nbatch(bdat))) {
    dat <- bigdat_batch(bdat, i)
    
    expect_equal(nrow(dat), nrows)
  }
})

test_that("basic example (bigmatrix)", {
  stopifnot(require(bigmemory))
  
  nrows <- 10
  ncols <- 5
  
  bmat <- list(as.big.matrix(matrix(1, nrows, ncols)), 
    as.big.matrix(matrix(2, nrows, ncols)))
  bdat <- bigdat(bmat, batch_size = 2)
  
  for(i in seq(1, bigdat_nbatch(bdat))) {
    dat <- bigdat_batch(bdat, i)
    
    expect_equal(nrow(dat), nrows)
  }
})
