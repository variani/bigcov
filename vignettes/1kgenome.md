# Analysis of 1K Genome data
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

## Data paths

We define a path to 1K dataset based on the output of `Sys.info()[["nodename"]]`.


```r
nodename <- Sys.info()[["nodename"]]

dpath <- switch(nodename,
  "tau" = "~/Data/1KGenome/",
  "your_nodename" = "/udd/redpr/mixed_model/1kG/",
  stop("unknown nodename"))
```


# GRM on a subset of common variants


```r
f <- file.path(dpath, "phase3_common_EUR503_pruned_100ksubset_bed.bed")
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/home/andrey/Data/1KGenome//phase3_common_EUR503_pruned_100ksubset_bed.bed]
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
pop <- read.table(file.path(dpath, "pop_label_phase3_2504_without_rels"), 
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
 45.256   1.804  47.180 
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

![](1kgenome_files/figure-html/com-1.png)<!-- -->

# GRM on a subset of rare variants


```r
f <- file.path(dpath, "phase3_rare_EUR503_nosingle_100ksubset_bed.bed")
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/home/andrey/Data/1KGenome//phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
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
pop <- read.table(file.path(dpath, "pop_label_phase3_2504_without_rels"), 
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
 44.808   1.632  46.665 
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

![](1kgenome_files/figure-html/rare-1.png)<!-- -->

# Jacard on a subset of rare variants 


```r
f <- file.path(dpath, "phase3_rare_EUR503_nosingle_100ksubset_bed.bed")
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/home/andrey/Data/1KGenome//phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
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
pop <- read.table(file.path(dpath, "pop_label_phase3_2504_without_rels"), 
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
 11.616   0.668  12.492 
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

![](1kgenome_files/figure-html/rare_jacard-1.png)<!-- -->


# Jacard on a subset of rare variants (<0.5%)


```r
f <- file.path(dpath, "phase3_rare_EUR503_nosingle_100ksubset_bed.bed")
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/home/andrey/Data/1KGenome//phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
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
pop <- read.table(file.path(dpath, "pop_label_phase3_2504_without_rels"), 
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
 12.312   0.100  12.461 
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

![](1kgenome_files/figure-html/rare_jacard2-1.png)<!-- -->


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
[1] bigcov_0.1.1    dplyr_0.5.0     magrittr_1.5    RSpectra_0.12-0
[5] BEDMatrix_1.4.0 rmarkdown_1.3   knitr_1.15.1    devtools_1.12.0

loaded via a namespace (and not attached):
 [1] codetools_0.2-15      digest_0.6.10         htmltools_0.3.5      
 [4] R6_2.2.0              assertthat_0.1        rprojroot_1.1        
 [7] grid_3.3.3            stringr_1.1.0         bigmemory.sri_0.1.3  
[10] testthat_1.0.2        tibble_1.2            lattice_0.20-34      
[13] DBI_0.5-1             foreach_1.4.3         roxygen2_5.0.1       
[16] iterators_1.0.8       Matrix_1.2-7.1        plyr_1.8.4           
[19] crochet_1.0.0         data.table_1.10.5     stringi_1.1.2        
[22] evaluate_0.10         yaml_2.1.14           tools_3.3.3          
[25] parallel_3.3.3        withr_1.0.2           synchronicity_1.1.9.1
[28] lazyeval_0.2.0        crayon_1.3.2          backports_1.0.4      
[31] memoise_1.0.0         bigmemory_4.5.19      Rcpp_0.12.11         
[34] doParallel_1.0.10    
```

# License

This document is licensed under the Creative Commons Attribution 4.0 International Public License. 

[Creative Commons License](http://creativecommons.org/licenses/by/4.0/)

# References
