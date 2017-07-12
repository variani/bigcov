context("parallel")

test_that("basic example (BEDMatrix)", {
  stopifnot(require(parallel))
  stopifnot(require(BEDMatrix))
  
  cores <- detectCores()
  if(cores > 1) {
    path <- system.file("extdata", "example.bed", package = "BEDMatrix")
    bmat <- BEDMatrix(path)
  
    t1 <- system.time(grm1 <- bigdat_grm(bmat, num_batches = 2, cores = 1))
    t2 <- system.time(grm2 <- bigdat_grm(bmat, num_batches = 2, cores = 2))
    
    expect_lt(t1[["elapsed"]], t2[["elapsed"]])
  }
})

