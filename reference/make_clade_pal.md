# Build a named Brewer palette with optional cycling behavior

Creates up to `n` visually distinct colors by sampling (and, when
needed, interpolating) from an RColorBrewer palette, optionally
reordering each block in a zig-zag pattern and progressively lightening
subsequent cycles.

## Usage

``` r
make_clade_pal(
  n,
  labels = NULL,
  pal = "Set3",
  alternating = FALSE,
  cycle_len = 0L,
  cycle_shift = 0.15
)
```

## Arguments

- n:

  Integer number of colors to generate.

- labels:

  Optional character vector of names to assign to the colors. If `NULL`,
  sequential integers are used.

- pal:

  RColorBrewer palette name provided to
  [`RColorBrewer::brewer.pal`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).

- alternating:

  Logical; when `TRUE`, reorders each chunk via a zig-zag pattern to
  maximize separation of adjacent hues.

- cycle_len:

  Positive integer forcing palette generation in repeating chunks of
  this size; `0` (default) disables chunking.

- cycle_shift:

  Scalar in `[0, 1]` controlling how strongly later chunks blend toward
  white to maintain distinguishability.

## Value

Named character vector of hex colors.
