# Optimize tree topology using C++ NNI moves

Performs nearest-neighbor interchange (NNI) hill-climbing to find the
tree topology that maximizes the belief-propagation score. Uses compiled
C++ routines for speed.

## Usage

``` r
optimize_tree_cpp(
  tree_init = NULL,
  logP,
  logA,
  max_iter = 100,
  outfile = NULL,
  resume = FALSE,
  ncores = 1,
  trace_interval = 5
)
```

## Arguments

- tree_init:

  A rooted `phylo` object used as the starting tree. Ignored when
  resuming from an existing trace.

- logP:

  A list of log-probability vectors (one per locus), as returned by
  [`convert_logliks_to_logP_list()`](https://teng-gao.github.io/mitodrift/reference/convert_logliks_to_logP_list.md)
  or `convert_logliks_to_logP_list_colmajor()`.

- logA:

  A numeric vector (or list of vectors) of log transition probabilities,
  flattened column-major from the transition matrix.

- max_iter:

  Integer; maximum number of NNI iterations.

- outfile:

  Optional file path for saving the tree trace (qs2 format).

- resume:

  Logical; if `TRUE` and `outfile` exists, resume from the last saved
  tree instead of starting fresh.

- ncores:

  Integer; number of threads for parallel NNI scoring.

- trace_interval:

  Integer; save the trace to `outfile` every this many iterations.

## Value

A `multiPhylo` list of trees visited during optimization, each carrying
a `logZ` element with the log-partition-function score.
