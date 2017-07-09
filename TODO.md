## TODO

- `bigdat_grm` function
    - create a shared memory matrix (bigmemory) to store/update results from each iteration/batch
    - merge `bigdat_jacard` into `bigdat_grm`
    - add more info. into computation output
        - the number of filtered SNPs: monomorphic, maf filter, NA
    - add `maf_max` argument, e.g. jacard for MAF bin [0.001; 0.005]
