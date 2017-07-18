context("jacard")

test_that("basic example: 0/1 entries", {
  nrows <- 5
  ncols <- 100
  
  # sparse matrix
  set.seed(1)
  vals <- rbinom(100, 1, 0.1)
  mat <- matrix(vals, nrows, ncols)
  
  jac1 <- jacard(mat)
  jac2 <- bigdat_grm(mat, "Jacard", filter = FALSE, num_batches = 2)
  
  expect_equal(as.numeric(jac1), as.numeric(jac2))
  expect_true(all(diag(jac1) == 1))
  expect_true(all(diag(jac2) == 1))

  # more dense matrix
  set.seed(1)
  vals <- rbinom(100, 1, 0.5)
  mat <- matrix(vals, nrows, ncols)
  
  jac3 <- jacard(mat)
  jac4 <- bigdat_grm(mat, "Jacard", filter = FALSE, num_batches = 2)
  
  expect_equal(as.numeric(jac3), as.numeric(jac4))
  expect_true(all(diag(jac3) == 1))
  expect_true(all(diag(jac4) == 1))
})

test_that("basic example: 0/1/2 entries", {
  nrows <- 5
  ncols <- 100
  
  # sparse matrix
  set.seed(1)
  vals1 <- rbinom(nrows*ncols, 1, 0.1)
  mat1 <- matrix(vals1, nrows, ncols)
  
  jac1 <- jacard(mat1)

  vals2 <- vals1
  vals2[vals2 == 1] <- 2
  mat2 <- matrix(vals2, nrows, ncols)

  jac2 <- jacard(mat2)

  expect_equal(as.numeric(jac1), as.numeric(jac2))
  expect_true(all(diag(jac1) == 1))
  expect_true(all(diag(jac2) == 1))

})


