# Moving-average smoothing across matrix rows

Applies a fixed-width moving average independently to every row of a
matrix. Missing values can be optionally ignored, and the `edge`
argument controls how windows at the matrix boundaries are handled. With
`edge = "partial"` (default), windows that extend past the matrix edge
shrink to the available cells, while `edge = "full"` enforces the full
window size by returning `NA` when a complete window cannot be formed.

## Usage

``` r
row_smooth(M, k, na_rm = FALSE, edge = c("partial", "full"))
```

## Arguments

- M:

  Numeric matrix to smooth (rows are independent series).

- k:

  Integer window width (must be \>= 1).

- na_rm:

  Logical; if `TRUE`, compute means ignoring `NA`s, otherwise only
  windows containing no `NA` values contribute.

- edge:

  Either `"partial"` or `"full"`; governs boundary handling as described
  above.

## Value

Matrix of the same dimensions as `M` containing smoothed values.
