# Decode a tree using CRF belief propagation (R version)

Constructs a conditional random field (CRF) on the tree with a single
shared transition matrix and computes per-variant marginal beliefs via
tree belief propagation. Optionally returns posterior means, MAP
assignments, or the full CRF objects.

## Usage

``` r
decode_tree(
  tn,
  A,
  liks,
  post_max = FALSE,
  store_bels = FALSE,
  store_crfs = FALSE,
  debug = FALSE,
  score_only = FALSE
)
```

## Arguments

- tn:

  A `phylo` object representing the tree topology.

- A:

  Transition matrix (square, `k x k`) with VAF bin midpoints as
  row/column names.

- liks:

  Named list of likelihood matrices (one per variant), each
  `k x n_cells`.

- post_max:

  Logical; if `TRUE`, also compute MAP (Viterbi) decoding.

- store_bels:

  Logical; if `TRUE`, store per-variant node and edge beliefs in the
  output.

- store_crfs:

  Logical; if `TRUE`, store a copy of the CRF object for each variant.

- debug:

  Logical; if `TRUE`, return a detailed list instead of just the
  `tbl_graph`.

- score_only:

  Logical; if `TRUE`, skip posterior computation and only attach
  log-partition scores.

## Value

A `tbl_graph` tree with per-variant posterior means (columns
`p_<variant>`) and a `logZ` vector of log-partition-function values.
When `debug = TRUE`, a list with additional diagnostic components.
