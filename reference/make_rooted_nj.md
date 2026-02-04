# Make a rooted NJ tree

Make a rooted NJ tree

## Usage

``` r
make_rooted_nj(vmat, dist_method = "manhattan", ncores = 1)
```

## Arguments

- vmat:

  A matrix of cell-by-variable values

- dist_method:

  The distance method to use

- ncores:

  Number of threads for
  [`parallelDist::parDist`](https://rdrr.io/pkg/parallelDist/man/parDist.html)
  (default: 1)

## Value

A phylo object
