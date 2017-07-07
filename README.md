# bigcov

[![travis-ci build status](https://travis-ci.org/variani/bigcov.svg?branch=master)](https://travis-ci.org/variani/bigcov)

Covariance/Correlation for big data

```
### timing #1
mat <- matrix(rnorm(1e6), nrow = 1e2, ncol = 1e4)

system.time(cov(mat)) # elapsed 6.834

system.time(bigcov(mat, num_batches = 2)) # elapsed 2.449

system.time(bigcov(mat, num_batches = 10)) # elapsed 1.351

### timing #2
mat <- matrix(rnorm(1e7), nrow = 1e5, ncol = 1e2)

system.time(cov(mat)) # elapsed 0.625

system.time(bigcov(mat, num_batches = 2)) # 0.372
```

## Data formats

- matrix
- bigmemory
- https://github.com/QuantGen/BEDMatrix (to be implemented)
