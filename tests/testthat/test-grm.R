context("grm")

test_that("basic example (BEDMatrix)", {
  stopifnot(require(BEDMatrix))
  
  path <- system.file("extdata", "example.bed", package = "BEDMatrix")
  bmat <- BEDMatrix(path)
  
  grm1 <- bigdat_grm(bmat, num_batches = 2, check_na = TRUE)
  grm2 <- bigdat_grm(bmat, num_batches = 2, check_na = FALSE)
})
