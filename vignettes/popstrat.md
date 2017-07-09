# Population stratification using GRM and Jacard
Andrey Ziyatdinov  
`r Sys.Date()`  



# About

Population stratification is an important part of association studies
conducted in a population consisting of different sub-populations.
The standard approach exploited in the analysis of common genetic markers is 
to estimate generic relatedness matrix (GRM) [@Patterson2006].
Recently, [@Prokopenko2015] proposed to use Jacard similarity matrix 
in the analysis of rare variants.
Here, we show practical examples on population stratification using R package [bigcov](https://github.com/variani/bigcov)
in the following settings.

| Simulation | MAF range | Number of samples | Number of Markers | Relationship |
|--|--|--|--|--|
| Common SNPs | [0.05; 0.5] | 200 | 1,000 | GRM |
| Rare SNPs | [0.001; 0.005] | 200 | 1,000 | GRM | 
| Rare SNPs | [0.001; 0.005] | 200 | 1,000, 2,000 & 4,000 | Jacard | 

We will use the same toy example of two population as in R package [jacpop](https://cran.r-project.org/web/packages/jacpop) 
reference manual (two populations with Fst = 0.1). 
Also, we will show an example of analysis for data stored in plink format.

## Work-horse functions in `bigcov`

Functions `bigdat_grm` and `bigdat_jacard` take input genotype data and 
convert them in `bigdat` format (also part of `bigcov`), which has a few practical features:

- support of several formats, including `matrix`, 
  `big.matrix` from R package [bigmemory](https://cran.r-project.org/web/packages/bigmemory/)
  and plink via R package [BEDMatrix](https://github.com/QuantGen/BEDMatrix);
- avoid loading the whole data matrix into RAM (for the later two formats) 
  and split the data into batches (either `num_batches` or `batch_size` arguments);
      - splitting into batches improves the performance, as shown by [SNPRelate](http://corearray.sourceforge.net/tutorials/SNPRelate/)
- support of several `bigmemory` or plink files, e.g. files per chromosomes;
- parallel computing (coming soon).
  
# Preliminaries


```r
library(BEDMatrix) # to read plink data efficiently

library(RSpectra) # to perform PCA with a few components, e.g. 2

#library(bigcov)
library(devtools)
load_all("~/git/variani/bigcov/")
```

## Function to simulate two populations

A toy example is adopted from https://cran.r-project.org/web/packages/jacpop.


```r
sim_pop <- function(N = 200, M = 1000, Fst = 0.1, maf_max = 0.5, maf_min = 0.05, 
  seed = 1)
{
  set.seed(seed)

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
  
  return(gdat)
}
```

## Function to plot PCA


```r
plot_pop <- function(mod, labs, ...)
{
  if(missing(labs)) {
    labs <- factor(rep("Pop", nrow(mod$vectors)))
  }
  
  plot(mod$vectors[, 1], mod$vectors[, 2], type = "n", 
    xlab = "PC1", ylab = "PC2", ...)
  text(mod$vectors[, 1], mod$vectors[, 2], label = labs, col = as.numeric(labs))
}
```

# Default simulation parameters


```r
N <- 200 # the number of individuals per population
M <- 1e3 # the number of SNPs
Fst <- 0.01
```

# Common SNPs


```r
# simulated genotype data
gdat <- sim_pop(N = N, M = M, Fst = Fst, maf_max = 0.5, maf_min = 0.05, seed = 1)
  
# compute GRM
A <- bigdat_grm(gdat, check_na = FALSE, num_batches = 2, verbose = 2)
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 2 batches
 - batch 1 / 2 
  -- # NAs 0 
 - batch 2 / 2 
  -- # NAs 0 
```

```r
# copmute PCA on GRM
mod <- eigs(A, k = 2)

# plot
labs <- factor(c(rep("Pop 1", N), rep("Pop 2", N)))
plot_pop(mod, labs, main = "GRM on common SNPs")
```

![](popstrat_files/figure-html/grm-1.png)<!-- -->

# Rare SNPs

## Population stratification using GRM


```r
# simulated genotype data
gdat <- sim_pop(N = N, M = M, Fst = Fst, maf_max = 0.005, maf_min = 0.001, seed = 1)
  
# compute GRM
A <- bigdat_grm(gdat, check_na = FALSE, num_batches = 2, verbose = 2)
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 2 batches
 - batch 1 / 2 
  -- filtered mono. markers: 206 
  -- # NAs 0 
 - batch 2 / 2 
  -- filtered mono. markers: 211 
  -- # NAs 0 
```

```r
# copmute PCA on GRM
mod <- eigs(A, k = 2)

# plot
labs <- factor(c(rep("Pop 1", N), rep("Pop 2", N)))
plot_pop(mod, labs, main = "GRM on rare SNPs")
```

![](popstrat_files/figure-html/grm_rare-1.png)<!-- -->

## Population stratification using Jacard


```r
# simulated genotype data
# - the same data as in the previous example, as the seed value is the same (1)
gdat <- sim_pop(N = N, M = M, Fst = Fst, maf_max = 0.005, maf_min = 0.001, seed = 1)

# compute Jacard
A <- bigdat_jacard(gdat, num_batches = 2, verbose = 2)
```

```
 - bigdat_jacard: computing: 2 batches
 - batch 1 / 2 
 - batch 2 / 2 
```

```r
# copmute PCA on GRM
mod <- eigs(A, k = 2)

# plot
labs <- factor(c(rep("Pop 1", N), rep("Pop 2", N)))
plot_pop(mod, labs, main = "Jacard on rare SNPs")
```

![](popstrat_files/figure-html/jacard_rare-1.png)<!-- -->

## Population stratification using Jacard and more SNPs


```r
# simulated genotype data
# - the simulated data is different, as we double the number of SNPs
M2 <- 2*M
gdat <- sim_pop(N = N, M = M2, Fst = Fst, maf_max = 0.005, maf_min = 0.001, seed = 1)

# compute Jacard
A <- bigdat_jacard(gdat, num_batches = 2, verbose = 2)
```

```
 - bigdat_jacard: computing: 2 batches
 - batch 1 / 2 
 - batch 2 / 2 
```

```r
# copmute PCA on GRM
mod <- eigs(A, k = 2)

# plot
labs <- factor(c(rep("Pop 1", N), rep("Pop 2", N)))
plot_pop(mod, labs, main = "Jacard on rare SNPs")
```

![](popstrat_files/figure-html/jacard_rare_M2-1.png)<!-- -->

## Population stratification using Jacard and even more SNPs


```r
M4 <- 4*M
gdat <- sim_pop(N = N, M = M4, Fst = Fst, maf_max = 0.005, maf_min = 0.001, seed = 1)

# compute Jacard
A <- bigdat_jacard(gdat, num_batches = 2, verbose = 2)
```

```
 - bigdat_jacard: computing: 2 batches
 - batch 1 / 2 
 - batch 2 / 2 
```

```r
# copmute PCA on GRM
mod <- eigs(A, k = 2)

# plot
labs <- factor(c(rep("Pop 1", N), rep("Pop 2", N)))
plot_pop(mod, labs, main = "Jacard on rare SNPs")
```

![](popstrat_files/figure-html/jacard_rare_M4-1.png)<!-- -->

# Example with plink files

We use dummy data in plink format from R package [BEDMatrix](https://github.com/QuantGen/BEDMatrix).


```r
path <- system.file("extdata", "example.bed", package = "BEDMatrix")
bmat <- BEDMatrix(path)

bmat
```

```
BEDMatrix: 50 x 1000 [/home/andrey/R/i686-pc-linux-gnu-library/3.3/BEDMatrix/extdata/example.bed]
```

```r
dim(bmat)
```

```
[1]   50 1000
```

```r
bmat[1:5, 1:5]
```

```
          snp0_A snp1_C snp2_G snp3_G snp4_G
per0_per0      0      1      1      1      0
per1_per1      1      1      1      1     NA
per2_per2      1      0      0      2      0
per3_per3      2      0      0      0      1
per4_per4      0      1      0      0      0
```


```r
N <- nrow(bmat)

grm <- bigdat_grm(bmat, num_batches = 2, verbose = 2)
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 2 batches
 - batch 1 / 2 
  -- # NAs 0 
 - batch 2 / 2 
  -- # NAs 0 
```

```r
mod <- eigs(grm, k = 2)

labs <- factor(rep("Dummy Pop", N))
plot_pop(mod, labs, main = "Dummy Data in plink format")
```

![](popstrat_files/figure-html/grm_plink-1.png)<!-- -->

# Future work

We need to test examples not covered in this document:

- real genotype data, e.g. 100K individual and 10M genotypes;
- parallel computing and splitting into many batches.

That will give us an idea whether R is powerful enough for big data in genetics
(spoiler: yes, it is). 

# R session info


```r
sessionInfo()
```

```
R version 3.3.3 (2017-03-06)
Platform: i686-pc-linux-gnu (32-bit)
Running under: Ubuntu 14.04.5 LTS

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] bigcov_0.1.1    readr_1.0.0     dplyr_0.5.0     magrittr_1.5   
[5] RSpectra_0.12-0 BEDMatrix_1.4.0 rmarkdown_1.3   knitr_1.15.1   
[9] devtools_1.12.0

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.8         bigmemory.sri_0.1.3 roxygen2_5.0.1     
 [4] lattice_0.20-34     R6_2.2.0            bigmemory_4.5.19   
 [7] stringr_1.1.0       plyr_1.8.4          tcltk_3.3.3        
[10] tools_3.3.3         grid_3.3.3          data.table_1.10.0  
[13] DBI_0.5-1           withr_1.0.2         htmltools_0.3.5    
[16] lazyeval_0.2.0      yaml_2.1.14         rprojroot_1.1      
[19] digest_0.6.10       assertthat_0.1      tibble_1.2         
[22] crayon_1.3.2        Matrix_1.2-7.1      codetools_0.2-15   
[25] testthat_1.0.2      memoise_1.0.0       evaluate_0.10      
[28] stringi_1.1.2       backports_1.0.4     crochet_1.0.0      
```

# License

This document is licensed under the Creative Commons Attribution 4.0 International Public License. 

[![Creative Commons License](http://i.creativecommons.org/l/by/4.0/88x31.png)](http://creativecommons.org/licenses/by/4.0/)

# References
