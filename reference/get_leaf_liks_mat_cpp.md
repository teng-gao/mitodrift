# Leaf likelihood matrix via C++ backend

Leaf likelihood matrix via C++ backend

## Usage

``` r
get_leaf_liks_mat_cpp(amat, dmat, vafs, eps = 0, ncores = 1L, log = FALSE)
```

## Arguments

- amat:

  Integer matrix of alt counts (variants x cells).

- dmat:

  Integer matrix of depths (same dimensions as `amat`).

- vafs:

  Numeric vector of VAF grid points (names optional).

- eps:

  Variant detection error rate to add to each VAF bin (default 0).

- ncores:

  Number of threads to use (default 1).

- log:

  Whether to return log-likelihoods instead of probabilities.

## Value

List of matrices, one per variant, with rows = VAF bins and columns =
cells.
