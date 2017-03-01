### data
dat <- iris[, -5]

cov(dat)

bigcov(dat)

### timing
mat <- matrix(rnorm(1e6), nrow = 1e2, ncol = 1e4)

system.time(cov(mat)) # elapsed 6.834

system.time(bigcov(mat, num_batches = 2)) # elapsed 2.449

system.time(bigcov(mat, num_batches = 10)) # elapsed 1.351

### timing #2
mat <- matrix(rnorm(1e7), nrow = 1e5, ncol = 1e2)

system.time(cov(mat)) # elapsed 0.625

system.time(bigcov(mat, num_batches = 2)) # 0.372
