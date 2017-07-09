context("jacard")

test_that("basic example", {
  nrows <- 5
  ncols <- 100
  
  # sparse matrix
  set.seed(1)
  vals <- rbinom(100, 1, 0.1)
  mat <- matrix(vals, nrows, ncols)
  
  jac1 <- jacard(mat)
  jac2 <- bigdat_jacard(mat, num_batches = 2)
  
  expect_equal(as.numeric(jac1), as.numeric(jac2))

  # more dense matrix
  set.seed(1)
  vals <- rbinom(100, 1, 0.5)
  mat <- matrix(vals, nrows, ncols)
  
  jac3 <- jacard(mat)
  jac4 <- bigdat_jacard(mat, num_batches = 2)
  
  expect_equal(as.numeric(jac3), as.numeric(jac4))
})


