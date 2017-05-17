context("tcrossprod")

test_that("tcrossprod is symmetric", {
  library(pls)
  
  data(gasoline)
  mat <- t(gasoline$NIR)
  
  cmat <- big_tcrossprod(mat)
  
  expect_true(Matrix::isSymmetric(cmat))
})

test_that("crossprod vs big_crossprod functions (4x4)", {
  data(iris)
  dat <- t(iris[, -5])
  
  mat <- tcrossprod(as.matrix(dat))
  b_mat <- big_tcrossprod(dat)
  
  expect_equal(as.numeric(cov_mat), as.numeric(bigcov_mat), tol = 1e-10)
})

test_that("crossprod vs big_crossprod functions (150x150)", {
  data(iris)
  dat <- iris[, -5]
  
  mat <- tcrossprod(as.matrix(dat))
  b_mat <- big_tcrossprod(dat)
  
  expect_equal(as.numeric(cov_mat), as.numeric(bigcov_mat), tol = 1e-10)
})
