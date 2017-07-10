context("grm")

test_that("basic example (BEDMatrix)", {
  stopifnot(require(BEDMatrix))
  
  path <- system.file("extdata", "example.bed", package = "BEDMatrix")
  bmat <- BEDMatrix(path)
  
  grm1 <- bigdat_grm(bmat, num_batches = 2, check_na = TRUE)
  grm2 <- bigdat_grm(bmat, num_batches = 2, check_na = FALSE)
})

test_that("all markers filtered out", {
  nrows <- 5
  ncols <- 100
  
  # sparse matrix
  set.seed(1)
  vals <- rbinom(100, 1, 0.5)
  mat <- matrix(vals, nrows, ncols)
  
  grm <- bigdat_grm(mat, num_batches = 2, maf_min = 0.6)
})

test_that("grm: jacard", {
  nrows <- 5
  ncols <- 100
  
  # sparse matrix
  set.seed(1)
  vals <- rbinom(100, 1, 0.1)
  mat <- matrix(vals, nrows, ncols)
  
  jac1 <- jacard(mat)
  jac2 <- bigdat_grm(mat, "Jacard", num_batches = 2)
  
  expect_equal(as.numeric(jac1), as.numeric(jac2))

  # more dense matrix
  set.seed(1)
  vals <- rbinom(100, 1, 0.5)
  mat <- matrix(vals, nrows, ncols)
  
  jac3 <- jacard(mat)
  jac4 <- bigdat_grm(mat, "Jacard", num_batches = 2)
  
  expect_equal(as.numeric(jac3), as.numeric(jac4))
})
