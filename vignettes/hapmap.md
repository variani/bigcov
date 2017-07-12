# Analysis of HapMap data
Andrey Ziyatdinov  
`r Sys.Date()`  


  
# Preliminaries


```r
library(utils) # for downloads
library(BEDMatrix) # to read plink data efficiently

library(RSpectra) # to perform PCA with a few components, e.g. 2

library(magrittr) 
library(dplyr) 

#library(bigcov)
library(devtools)
load_all("~/git/variani/bigcov/")
```


## Function to plot PCA


```r
plot_pop <- function(mod, labs, ...)
{
  plot(mod$vectors[, 1], mod$vectors[, 2], type = "n", 
    xlab = "PC1", ylab = "PC2", ...)
  text(mod$vectors[, 1], mod$vectors[, 2], label = labs, col = as.numeric(labs))
}
```

## Download and extract data 

http://zzz.bwh.harvard.edu/plink/res.shtml#hapmap

>The HapMap genotype data (the latest is release 23) are available here as PLINK binary filesets. The SNPs are currently coded according NCBI build 36 coordinates on the forward strand. Several versions are available here: the entire dataset (a single, very large fileset: you will need a computer with at least 2Gb of RAM to load this file).
>
>The filtered SNP set refers to a list of SNPs that have MAF greater than 0.01 and genotyping rate greater than 0.95 in the 60 CEU founders. This fileset is probably a good starting place for imputation in samples of European descent. Filtered versions of the other HapMap panels will be made available shortly.


```r
dir_data <- "data-hapmap" 
dir.create(dir_data, showWarnings = FALSE)

# description file
url_pop <- "http://zzz.bwh.harvard.edu/plink/dist/hapmap.pop"
fn_pop <- basename(url_pop)
download.file(url_pop, file.path(dir_data, fn_pop))

# genotype file
url_gen <- "http://zzz.bwh.harvard.edu/plink/dist/hapmap_JPT_CHB_r23a.zip"
fn_gen <- basename(url_gen)
download.file(url_gen, file.path(dir_data, fn_gen))
```


```r
unzip(file.path(dir_data, fn_gen), exdir = dir_data)
```


```r
list.files(dir_data)
```

```
[1] "hapmap_JPT_CHB_r23a.bed" "hapmap_JPT_CHB_r23a.bim"
[3] "hapmap_JPT_CHB_r23a.fam" "hapmap_JPT_CHB_r23a.zip"
[5] "hapmap.pop"             
```


## Read files by BEDMatrix


```r
bed_gen <- paste0(strsplit(fn_gen, '[.]')[[1]][1], ".bed")
system.time({
  bmat <- BEDMatrix(file.path(dir_data, bed_gen))
})
```

```
   user  system elapsed 
  3.683   0.088   4.989 
```


```r
bmat
```

```
BEDMatrix: 90 x 3998895 [data-hapmap/hapmap_JPT_CHB_r23a.bed]
```

```r
dim(bmat)
```

```
[1]      90 3998895
```

### Extract IDs of individuals


```r
ids <- bmat@dnames[[1]]
head(ids)
```

```
[1] "NA18524_NA18524" "NA18526_NA18526" "NA18529_NA18529" "NA18532_NA18532"
[5] "NA18537_NA18537" "NA18540_NA18540"
```

```r
ids <- strsplit(ids, "_") %>% sapply(function(x) x[1])
head(ids)
```

```
[1] "NA18524" "NA18526" "NA18529" "NA18532" "NA18537" "NA18540"
```


```r
pop <- read.table(file.path(dir_data, fn_pop), header = FALSE)
names(pop) <- c("name", "id", "pop")

pop <- left_join(data_frame(id = ids), pop, by = "id")
stopifnot(nrow(pop) == length(ids))
```

## Summary on data set 

| Characteristic | Value |
|--|--|
| File | `hapmap_JPT_CHB_r23a.zip` |
| Description | [JPT+CHB (release 23, 90 individuals, 3.99 million SNPs)](http://zzz.bwh.harvard.edu/plink/res.shtml#hapmap) |
| Orignial file | `hapmap_JPT_CHB_r23a.bed` |
| File format in R | `BEDMatrix` |
| Number of individuals | 90 |
| Number of markers | 3,998,895 |

# GRM on all markers


```r
system.time({
  grm <- bigdat_grm(bmat, batch_size = 1e5, verbose = 1)
})
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 40 batches
 - clean markers used in the analysis: 2573397 / 3998895 
```

```
   user  system elapsed 
 39.069   2.817  42.029 
```


```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on all markers")
```

![](hapmap_files/figure-html/pca_all-1.png)<!-- -->

# GRM on common markers (> 5%)


```r
system.time({
  grm <- bigdat_grm(bmat, batch_size = 1e5, maf_min = 0.05, verbose = 1)
})
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 40 batches
 - clean markers used in the analysis: 2133712 / 3998895 
```

```
   user  system elapsed 
 27.435   2.696  29.666 
```


```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on common markers")
```

![](hapmap_files/figure-html/pca_common-1.png)<!-- -->

# Jacard on rare markers (<5%)

The genotype data does not have enough rare variants, e.g. < 1%.
Thus, the Jacard matrix is the identity matrix
if the analysis is run on genotypes with MAF below 1%.


```r
system.time({
  grm <- bigdat_grm(bmat, "Jacard", batch_size = 1e5, maf_max = 0.05, verbose = 1) 
})
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 40 batches
 - clean markers used in the analysis: 454802 / 3998895 
```

```
   user  system elapsed 
 13.057   1.494  14.561 
```



```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "Jacard on rare markers")
```

![](hapmap_files/figure-html/pca_jacard_rare-1.png)<!-- -->


# GRM on rare markers (<5%)


```r
system.time({
  grm <- bigdat_grm(bmat, batch_size = 1e5, maf_max = 0.05, verbose = 1) 
})
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 40 batches
 - clean markers used in the analysis: 454802 / 3998895 
```

```
   user  system elapsed 
 12.984   1.533  14.413 
```



```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on rare markers")
```

![](hapmap_files/figure-html/pca_grm_rare-1.png)<!-- -->

# Parallel computing 

We compute the previous example on a single or two CPUs.
The following code does not work on Mac with BLAS supporting multi threading.


```r
system.time(bigdat_grm(bmat, batch_size = 1e5, maf_max = 0.05))

system.time(bigdat_grm(bmat, batch_size = 1e5, maf_max = 0.05, cores = 2))
```
  
# R session info


```r
sessionInfo()
```

```
R version 3.4.0 (2017-04-21)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X El Capitan 10.11.6

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] bigcov_0.1.1    dplyr_0.5.0     magrittr_1.5    RSpectra_0.12-0
[5] BEDMatrix_1.4.0 rmarkdown_1.5   knitr_1.15.1    devtools_1.13.1

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.10          compiler_3.4.0        plyr_1.8.4           
 [4] iterators_1.0.8       tools_3.4.0           crochet_1.0.0        
 [7] testthat_1.0.2        digest_0.6.12         evaluate_0.10        
[10] memoise_1.1.0         tibble_1.3.0          lattice_0.20-35      
[13] Matrix_1.2-9          foreach_1.4.3         DBI_0.6-1            
[16] synchronicity_1.1.9.1 commonmark_1.2        yaml_2.1.14          
[19] parallel_3.4.0        withr_1.0.2           stringr_1.2.0        
[22] roxygen2_6.0.1        xml2_1.1.1            rprojroot_1.2        
[25] grid_3.4.0            data.table_1.10.4     R6_2.2.1             
[28] bigmemory_4.5.19      bigmemory.sri_0.1.3   backports_1.0.5      
[31] codetools_0.2-15      htmltools_0.3.6       assertthat_0.2.0     
[34] stringi_1.1.5         lazyeval_0.2.0        doParallel_1.0.10    
[37] crayon_1.3.2         
```

# License

This document is licensed under the Creative Commons Attribution 4.0 International Public License. 

[![Creative Commons License](http://i.creativecommons.org/l/by/4.0/88x31.png)](http://creativecommons.org/licenses/by/4.0/)

# References
