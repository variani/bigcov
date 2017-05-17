context("crossprod")

test_that("crossprod is symmetric", {
  library(pls)
  
  data(gasoline)
  mat <- gasoline$NIR
  
  cmat <- big_crossprod(mat)
  
  expect_true(Matrix::isSymmetric(cmat))
})

test_that("crossprod vs big_crossprod functions", {
  data(iris)
  dat <- iris[, -5]
  
  mat <- crossprod(as.matrix(dat))
  b_mat <- big_crossprod(dat)
  
  expect_equal(as.numeric(cov_mat), as.numeric(bigcov_mat), tol = 1e-10)
})

