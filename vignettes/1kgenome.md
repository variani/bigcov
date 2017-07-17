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

# GRM on a subset of common variants


```r
f <- "~/Data/1KGenome/phase3_common_EUR503_pruned_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/Users/andreyziyatdinov/Data/1KGenome/phase3_common_EUR503_pruned_100ksubset_bed.bed]
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
pop <- read.table("~/Data/1KGenome/pop_label_phase3_2504_without_rels", 
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
  6.007   0.661   6.607 
```

```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on a pre-selected subset of common markers")
```

![](1kgenome_files/figure-html/com-1.png)<!-- -->

# GRM on a subset of rare variants


```r
f <- "~/Data/1KGenome/phase3_rare_EUR503_nosingle_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/Users/andreyziyatdinov/Data/1KGenome/phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
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
pop <- read.table("~/Data/1KGenome/pop_label_phase3_2504_without_rels", 
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
  5.271   0.605   5.552 
```

```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on a pre-selected subset of rare markers")
```

![](1kgenome_files/figure-html/rare-1.png)<!-- -->

# Jacard on a subset of rare variants 


```r
f <- "~/Data/1KGenome/phase3_rare_EUR503_nosingle_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/Users/andreyziyatdinov/Data/1KGenome/phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
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
pop <- read.table("~/Data/1KGenome/pop_label_phase3_2504_without_rels", 
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
  4.132   0.958   5.103 
```

```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on a pre-selected subset of rare markers")
```

![](1kgenome_files/figure-html/rare_jacard-1.png)<!-- -->


# Jacard on a subset of rare variants (<0.5%)


```r
f <- "~/Data/1KGenome/phase3_rare_EUR503_nosingle_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

bmat
```

```
BEDMatrix: 503 x 100000 [/Users/andreyziyatdinov/Data/1KGenome/phase3_rare_EUR503_nosingle_100ksubset_bed.bed]
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
pop <- read.table("~/Data/1KGenome/pop_label_phase3_2504_without_rels", 
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
  4.789   0.766   5.825 
```

```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on a pre-selected subset of common markers (<0.05%)")
```

![](1kgenome_files/figure-html/rare_jacard2-1.png)<!-- -->

# GRM on a full set of common variants (>5%)


```r
bed <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bed"
fam <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.fam"
bim <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bim"

ret <- system(paste("wc -l", fam), intern = TRUE)
n <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

ret <- system(paste("wc -l", bim), intern = TRUE)
p <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

path <- path.expand(bed)
bmat <- BEDMatrix(path = path, p = p)

bmat
```

```
BEDMatrix: 503 x 24887325 [/Users/andreyziyatdinov/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bed]
```

```r
dim(bmat)
```

```
[1]      503 24887325
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
pop <- read.table("~/Data/1KGenome/pop_label_phase3_2504_without_rels", 
    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(pop = factor(pop)) %>% 
  left_join(data_frame(sample = ids), ., by = "sample")

system.time({
  grm <- bigdat_grm(bmat, batch_size = 5e4, maf_min = 0.05, verbose = 1)
})
```

```
  -- filters: callrate / mono / maf_min /  / check_na 
 - bigdat_tcrossprod: computing `tcrossprod`: 498 batches
 - clean markers used in the analysis: 7153139 / 24887325 
```

```
    user   system  elapsed 
 894.530  160.551 1062.520 
```

```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on a full set of common markers (>5%)")
```

![](1kgenome_files/figure-html/com_full-1.png)<!-- -->


# Jacard on a full set of rare variants (<1%)


```r
bed <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bed"
fam <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.fam"
bim <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bim"

ret <- system(paste("wc -l", fam), intern = TRUE)
n <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

ret <- system(paste("wc -l", bim), intern = TRUE)
p <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

path <- path.expand(bed)
bmat <- BEDMatrix(path = path, p = p)

bmat
```

```
BEDMatrix: 503 x 24887325 [/Users/andreyziyatdinov/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bed]
```

```r
dim(bmat)
```

```
[1]      503 24887325
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
pop <- read.table("~/Data/1KGenome/pop_label_phase3_2504_without_rels", 
    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(pop = factor(pop)) %>% 
  left_join(data_frame(sample = ids), ., by = "sample")

system.time({
  grm <- bigdat_grm(bmat, "Jacard", batch_size = 5e4, maf_max = 0.01, verbose = 1)
})
```

```
  -- filters: callrate / mono /  / maf_max / check_na 
 - bigdat_tcrossprod: computing `tcrossprod`: 498 batches
 - clean markers used in the analysis: 14813099 / 24887325 
```

```
    user   system  elapsed 
1049.317  215.654 1298.037 
```

```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "Jacard on a full set of rare markers (<1%)")
```

![](1kgenome_files/figure-html/rare_full-1.png)<!-- -->

# Jacard on a full set of rare variants (<0.5%)


```r
bed <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bed"
fam <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.fam"
bim <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bim"

ret <- system(paste("wc -l", fam), intern = TRUE)
n <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

ret <- system(paste("wc -l", bim), intern = TRUE)
p <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

path <- path.expand(bed)
bmat <- BEDMatrix(path = path, p = p)

bmat
```

```
BEDMatrix: 503 x 24887325 [/Users/andreyziyatdinov/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bed]
```

```r
dim(bmat)
```

```
[1]      503 24887325
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
pop <- read.table("~/Data/1KGenome/pop_label_phase3_2504_without_rels", 
    header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(pop = factor(pop)) %>% 
  left_join(data_frame(sample = ids), ., by = "sample")

system.time({
  grm <- bigdat_grm(bmat, "Jacard", batch_size = 1e5, maf_max = 0.005, verbose = 1)
})
```

```
  -- filters: callrate / mono /  / maf_max / check_na 
 - bigdat_tcrossprod: computing `tcrossprod`: 249 batches
 - clean markers used in the analysis: 13461033 / 24887325 
```

```
    user   system  elapsed 
1111.510  217.223 1364.575 
```

```r
grm[1:5, 1:5]
```

```
            1           2           3           4           5
1 1.164451233 0.002590149 0.002021923 0.002776386 0.003554834
2 0.002590149 1.003340013 0.001389468 0.007127107 0.005116180
3 0.002021923 0.001389468 1.003553779 0.012827612 0.007755358
4 0.002776386 0.007127107 0.012827612 1.006500702 0.013766705
5 0.003554834 0.005116180 0.007755358 0.013766705 1.097581595
```

```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "Jacard on a full set of rare markers (<0.1%)")
```

![](1kgenome_files/figure-html/rare_full2-1.png)<!-- -->

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

[Creative Commons License](http://creativecommons.org/licenses/by/4.0/)

# References
