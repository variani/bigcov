### inc
library(BEDMatrix)

library(magrittr)

### common
f <- "~/Data/1KGenome/phase3_common_EUR503_pruned_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

#grm <- bigdat_grm(bmat, batch_size = 1e4, verbose = 2)

### rare
f <- "~/Data/1KGenome/phase3_rare_EUR503_nosingle_100ksubset_bed.bed"
bmat <- BEDMatrix(f)

### all
bed <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bed"
fam <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.fam"
bim <- "~/Data/1KGenome/ALL_withDI_EUR503_phase3_nomultialleles_nodupls.bim"

ret <- system(paste("wc -l", fam), intern = TRUE)
n <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

ret <- system(paste("wc -l", bim), intern = TRUE)
p <- strsplit(ret, " ") %>% .[[1]] %>% as.integer %>% .[!is.na(.)]

path <- path.expand(bed)
bmat <- BEDMatrix(path = path, n = n, p = p)




