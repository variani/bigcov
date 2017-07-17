---
title: "Analysis of 1K Genome data"
author: "Andrey Ziyatdinov"
date: "2017-07-17"
output:
  html_document:
    theme: united
    toc: true
    keep_md: true
bibliography: ref.bib    
---


  
# Preliminaries


```r
library(utils) # for downloads
library(BEDMatrix) # to read plink data efficiently

library(RSpectra) # to perform PCA with a few components, e.g. 2

library(magrittr) 
library(dplyr) 

#library(bigcov)
library(devtools)
#load_all("~/git/variani/bigcov/")
load_all()
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

# GRM on a subset of common variants


```r
f <- "/udd/redpr/mixed_model/1kG/phase3_common_EUR503_pruned_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/udd/redpr/mixed_model/1kG/phase3_common_EUR503_pruned_100ksubset_bed.bed]
```

```r
dim(bmat)
```

```
[1]    503 100000
```

```r
ids <- rownames(bmat)
ids <- strsplit(ids, "_") %>% sapply(function(x) x[1])
head(ids)
```

```
[1] "HG00096" "HG00097" "HG00099" "HG00100" "HG00101" "HG00102"
```

```r
pop <- read.table("/udd/redpr/mixed_model/1kG/pop_label_phase3_2504_without_rels", 
    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(pop = factor(pop)) %>% 
  left_join(data_frame(sample = ids), ., by = "sample")

system.time({
  grm <- bigdat_grm(bmat, batch_size = 1e4, verbose = 1)
})
```

```
  -- filters: callrate / mono /  /  / check_na 
 - bigdat_tcrossprod: computing `tcrossprod`: 10 batches
 - clean markers used in the analysis: 100000 / 100000 
```

```
   user  system elapsed 
 44.721   0.287  45.068 
```

```r
mod <- eigs(grm, k = 2)
mod2<-prcomp(grm)
names(mod2)[2]<-'vectors'

labs <- pop$pop
par(mfrow=c(1, 2))
plot_pop(mod, labs, main = "GRM on common markers, using eigs")
plot_pop(mod2,labs,main='GRM on common markers, using base::prcomp')
```

![plot of chunk com](figure/com-1.png)

# GRM on a subset of rare variants


```r
f <- "/udd/redpr/mixed_model/1kG/phase3_rare_EUR503_nosingle_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/udd/redpr/mixed_model/1kG/phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
```

```r
dim(bmat)
```

```
[1]    503 100000
```

```r
ids <- rownames(bmat)
ids <- strsplit(ids, "_") %>% sapply(function(x) x[1])
head(ids)
```

```
[1] "HG00096" "HG00097" "HG00099" "HG00100" "HG00101" "HG00102"
```

```r
pop <- read.table("/udd/redpr/mixed_model/1kG/pop_label_phase3_2504_without_rels", 
    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(pop = factor(pop)) %>% 
  left_join(data_frame(sample = ids), ., by = "sample")

system.time({
  grm <- bigdat_grm(bmat, batch_size = 1e4, verbose = 1)
})
```

```
  -- filters: callrate / mono /  /  / check_na 
 - bigdat_tcrossprod: computing `tcrossprod`: 10 batches
 - clean markers used in the analysis: 100000 / 100000 
```

```
   user  system elapsed 
 48.138   0.214  48.439 
```

```r
mod <- eigs(grm, k = 2)
mod2<-prcomp(grm)
names(mod2)[2]<-'vectors'

labs <- pop$pop
par(mfrow=c(1, 2))
plot_pop(mod, labs, main = "GRM on rare markers, using eigs")
plot_pop(mod2,labs,main='GRM on rare markers, using base::prcomp')
```

![plot of chunk rare](figure/rare-1.png)

# Jacard on a subset of rare variants 


```r
f <- "/udd/redpr/mixed_model/1kG/phase3_rare_EUR503_nosingle_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/udd/redpr/mixed_model/1kG/phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
```

```r
dim(bmat)
```

```
[1]    503 100000
```

```r
ids <- rownames(bmat)
ids <- strsplit(ids, "_") %>% sapply(function(x) x[1])
head(ids)
```

```
[1] "HG00096" "HG00097" "HG00099" "HG00100" "HG00101" "HG00102"
```

```r
pop <- read.table("/udd/redpr/mixed_model/1kG/pop_label_phase3_2504_without_rels", 
    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(pop = factor(pop)) %>% 
  left_join(data_frame(sample = ids), ., by = "sample")

system.time({
  grm <- bigdat_grm(bmat, "Jacard", batch_size = 1e4, verbose = 1)
})
```

```
  -- filters: callrate / mono /  /  / check_na 
 - bigdat_tcrossprod: computing `tcrossprod`: 10 batches
 - clean markers used in the analysis: 100000 / 100000 
```

```
   user  system elapsed 
 14.287   0.217  14.531 
```

```r
mod <- eigs(grm, k = 2)
mod2<-prcomp(grm)
names(mod2)[2]<-'vectors'

labs <- pop$pop
par(mfrow=c(1, 2))
plot_pop(mod, labs, main = "JAC on rare markers, using eigs")
plot_pop(mod2,labs,main='JAC on rare markers, using base::prcomp')
```

![plot of chunk rare_jacard](figure/rare_jacard-1.png)


# Jacard on a subset of rare variants (<0.5%)


```r
f <- "/udd/redpr/mixed_model/1kG/phase3_rare_EUR503_nosingle_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/udd/redpr/mixed_model/1kG/phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
```

```r
dim(bmat)
```

```
[1]    503 100000
```

```r
ids <- rownames(bmat)
ids <- strsplit(ids, "_") %>% sapply(function(x) x[1])
head(ids)
```

```
[1] "HG00096" "HG00097" "HG00099" "HG00100" "HG00101" "HG00102"
```

```r
pop <- read.table("/udd/redpr/mixed_model/1kG/pop_label_phase3_2504_without_rels", 
    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(pop = factor(pop)) %>% 
  left_join(data_frame(sample = ids), ., by = "sample")

system.time({
  grm <- bigdat_grm(bmat, "Jacard", batch_size = 1e4, maf_max = 5e-3, verbose = 1)
})
```

```
  -- filters: callrate / mono /  / maf_max / check_na 
 - bigdat_tcrossprod: computing `tcrossprod`: 10 batches
 - clean markers used in the analysis: 75197 / 100000 
```

```
   user  system elapsed 
 14.807   0.051  14.884 
```

```r
mod <- eigs(grm, k = 2)
mod2<-prcomp(grm)
names(mod2)[2]<-'vectors'

labs <- pop$pop
par(mfrow=c(1, 2))
plot_pop(mod, labs, main = "JAC on rare markers, using eigs")
plot_pop(mod2,labs,main='JAC on rare markers, using base::prcomp')
```

![plot of chunk rare_jacard2](figure/rare_jacard2-1.png)


# R session info


```r
sessionInfo()
```

```
R version 3.4.0 (2017-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS release 6.8 (Final)

Matrix products: default
BLAS: /app/R-3.4.0_nonMKL@i86-rhel6.0/lib64/R/lib/libRblas.so
LAPACK: /app/R-3.4.0_nonMKL@i86-rhel6.0/lib64/R/lib/libRlapack.so

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
[1] bindrcpp_0.2    bigcov_0.1.1    devtools_1.13.2 dplyr_0.7.1    
[5] magrittr_1.5    RSpectra_0.12-0 BEDMatrix_1.4.0 markdown_0.8   
[9] knitr_1.16     

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.11          highr_0.6             compiler_3.4.0       
 [4] plyr_1.8.4            bindr_0.1             iterators_1.0.8      
 [7] tools_3.4.0           crochet_1.0.0         testthat_1.0.2       
[10] digest_0.6.12         evaluate_0.10.1       memoise_1.1.0        
[13] tibble_1.3.3          lattice_0.20-35       pkgconfig_2.0.1      
[16] rlang_0.1.1           Matrix_1.2-10         foreach_1.4.3        
[19] synchronicity_1.1.9.1 commonmark_1.2        parallel_3.4.0       
[22] withr_1.0.2           stringr_1.2.0         roxygen2_6.0.1       
[25] xml2_1.1.1            grid_3.4.0            data.table_1.10.4    
[28] glue_1.1.1            R6_2.2.2              bigmemory_4.5.19     
[31] bigmemory.sri_0.1.3   codetools_0.2-15      assertthat_0.2.0     
[34] stringi_1.1.5         doParallel_1.0.10     crayon_1.3.2         
```

# License

This document is licensed under the Creative Commons Attribution 4.0 International Public License. 

[Creative Commons License](http://creativecommons.org/licenses/by/4.0/)

# References
