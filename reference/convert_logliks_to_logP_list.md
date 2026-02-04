# Convert log-likelihood matrices to a log-probability list (row-major)

Convert log-likelihood matrices to a log-probability list (row-major)

## Usage

``` r
convert_logliks_to_logP_list(logliks, phy)
```

## Arguments

- logliks:

  Named list of log-likelihood matrices (one per variant), each of
  dimension `k x n_cells`.

- phy:

  A `phylo` object whose tip labels determine column ordering.

## Value

A named list of numeric vectors, each of length `k * n_nodes`.
