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

test_that("basic example (BEDMatrix)", {
  stopifnot(require(BEDMatrix))
  
  path <- system.file("extdata", "example.bed", package = "BEDMatrix")
  bmat <- BEDMatrix(path)
  nrows <- nrow(bmat)
  
  bdat <- bigdat(bmat, num_batches = 2)
    
  for(i in seq(1, bigdat_nbatch(bdat))) {
    dat <- bigdat_batch(bdat, i)
    
    expect_equal(nrow(dat), nrows)
  }
})

test_that("basic example (Gaston)", {
  stopifnot(require(gaston))
  
  path <- system.file("extdata", "example.bed", package = "BEDMatrix")
  bmat <- read.bed.matrix(path)
  nrows <- nrow(bmat)
  
  bdat <- bigdat(bmat, num_batches = 2)
    
  for(i in seq(1, bigdat_nbatch(bdat))) {
    dat <- bigdat_batch(bdat, i)
    
    expect_equal(nrow(dat), nrows)
  }
})

test_that("assoc example (Gaston)", {
  stopifnot(require(gaston))
  
  path <- system.file("extdata", "example.bed", package = "BEDMatrix")
  bmat <- read.bed.matrix(path)
  nrows <- nrow(bmat)
  
  bdat <- bigdat(bmat, batch_size = 10, as_matrix = FALSE)
  
  N <- bigdat_nrow1(bdat)
  y <- rnorm(N)

  # batch 1
  X <- bigdat_batch(bdat, 1)
  
  assoc <- gaston::association.test(X, y)

})
