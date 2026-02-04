# Generate VAF bin midpoints

Creates `k + 2` evenly spaced VAF bins spanning \[0, 1\] (including
boundary bins at 0 and 1) and returns their midpoints.

## Usage

``` r
get_vaf_bins(k)
```

## Arguments

- k:

  Integer; number of interior VAF bins. The total number of bins is
  `k + 2`.

## Value

Numeric vector of bin midpoints of length `k + 2`.
