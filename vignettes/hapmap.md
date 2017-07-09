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
  6.348   0.128   8.872 
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

# GRM on all markers


```r
system.time({
  grm <- bigdat_grm(bmat, batch_size = 1e5, verbose = 1)
})
```

```
 - bigdat_tcrossprod: computing `tcrossprod`: 40 batches
```

```
   user  system elapsed 
172.336   5.436 179.352 
```


```r
mod <- eigs(grm, k = 2)

labs <- pop$pop
plot_pop(mod, labs, main = "GRM on all markers")
```

![](hapmap_files/figure-html/pca_all-1.png)<!-- -->

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

![Creative Commons License](http://creativecommons.org/licenses/by/4.0/)

# References
