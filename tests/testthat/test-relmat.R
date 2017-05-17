context("relmat")

test_that("tcrossprod is symmetric", {
  center <- TRUE
  scale <- TRUE
  
  data(iris)
  dat <- t(iris[, -5])
  
  cov1 <- big_relmat(dat, center = center, scale = scale)
  cov2 <- big_cov(t(dat))

  cor1 <- cov2cor(cov1)
  cor2 <- cov2cor(cov2)
  
  diff <- abs(as.numeric(cor1) - as.numeric(cor2))
  diff <- mean(diff)
  
  expect_true(diff > 0.1)
})

