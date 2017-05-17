context("bigcov")

test_that("cov vs bigcov functions", {
  data(iris)
  dat <- iris[, -5]
  
  cov_mat <- cov(dat)
  bigcov_mat <- bigcov(dat)
  
  expect_equal(as.numeric(cov_mat), as.numeric(bigcov_mat), tol = 1e-10)
})

test_that("cov. is symmetric", {
  library(pls)
  
  data(gasoline)
  mat <- gasoline$NIR
  
  cmat <- bigcov(mat)
  
  expect_true(Matrix::isSymmetric(cmat))
})
