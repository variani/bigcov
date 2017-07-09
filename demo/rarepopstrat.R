# Toy example adopted from https://cran.r-project.org/web/packages/jacpop
# - simulate genotypes in 2 populations

### inc
library(RSpectra)

### par
N <- 2e2 # the number of individuals per population
M <- 1e3 # the number of SNPs
Fst <- 0.01

maf_max <- 0.005
maf_min <- 0.001

### simulate data
set.seed(1)

maf_values <- runif(M, maf_min, maf_max)

freq1 <- sapply(1:M, function(i) rbeta(1, 
  maf_values[i] * (1 - Fst) / Fst, 
  (1 - maf_values[i]) * (1 - Fst) / Fst))
freq2 <- sapply(1:M, function(i) rbeta(1, 
  maf_values[i] * (1 - Fst) / Fst, 
  (1 - maf_values[i]) * (1 - Fst) / Fst))
  
gdat1 <- sapply(1:M, function(i) sample(c(0, 1, 2), N, replace = TRUE,
  prob = c(((1 - freq1[i])^2), (2 * freq1[i] * (1 - freq1[i])), (freq1[i]^2))))
gdat2 <- sapply(1:M, function(i) sample(c(0, 1, 2), N, replace = TRUE,
  prob = c(((1 - freq2[i])^2), (2 * freq2[i] * (1 - freq2[i])), (freq2[i]^2))))

gdat <- rbind(gdat1, gdat2)

### compute GRM
A <- bigdat_grm(gdat, check_na = FALSE, verbose = 2)

### copmute PCA
mod <- eigs(A, k = 2)

### plot
labs <- factor(c(rep("Pop 1", N), rep("Pop 2", N)))

plot(mod$vectors[, 1], mod$vectors[, 2], type = "n", 
  xlab = "PC1", ylab = "PC2", main = "GRM on rare SNPs")
text(mod$vectors[, 1], mod$vectors[, 2], label = labs, col = as.numeric(labs))

### plot excluding an outlier (visually detected from the previous plot)
#ind <- which(mod$vectors[, 2] < 0.6)

#plot(mod$vectors[ind, 1], mod$vectors[ind, 2], type = "n", 
#  xlab = "PC1", ylab = "PC2", main = "GRM on rare SNPs (excluding outliers)")
#text(mod$vectors[ind, 1], mod$vectors[ind, 2], label = labs, col = as.numeric(labs))

