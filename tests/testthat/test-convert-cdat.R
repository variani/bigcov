context("convert cdat")

test_that("cov vs bigcov functions", {
  data(iris)
  dat <- iris[, -5]
  
  cmat <- cov(dat)
  
  cdat_all <- convert_cdat(cmat)
  cdat_upper_diag <- convert_cdat(cmat, lower = FALSE)
  cdat_upper <- convert_cdat(cmat, diag = FALSE, lower = FALSE)
  
  p <- ncol(cmat)
  
  expect_equal(nrow(cdat_all), p * p) # 16
  expect_equal(nrow(cdat_upper_diag), (p * p - p) / 2 + p) # 10
  expect_equal(nrow(cdat_upper), (p * p - p) / 2) # 6
})

