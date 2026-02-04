# Convert allele count matrix (amat) and total count matrix (dmat) to long format

Convert allele count matrix (amat) and total count matrix (dmat) to long
format

## Usage

``` r
mat_to_long(amat, dmat)
```

## Arguments

- amat:

  Allele count matrix with variants as rows and cells as columns

- dmat:

  Total count matrix with variants as rows and cells as columns

## Value

A data.table in long format with columns: variant, cell, a (allele
count), d (total count)
