# Run tree-topology MCMC in batches with convergence monitoring

Runs MCMC sampling in fixed-size batches, computing ASDSF convergence
diagnostics between batches. Supports automatic stopping when ASDSF
drops below a threshold and resume from a previous run.

## Usage

``` r
run_tree_mcmc_batch(
  phy_init,
  logP_list,
  logA_vec,
  outfile,
  diagfile = NULL,
  diag = TRUE,
  max_iter = 100,
  nchains = 1,
  ncores = 1,
  ncores_qs = 1,
  batch_size = 1000,
  conv_thres = NULL,
  resume = FALSE
)
```

## Arguments

- phy_init:

  A rooted `phylo` object used as the starting tree.

- logP_list:

  List of log-probability vectors (one per locus).

- logA_vec:

  Numeric vector of log transition probabilities.

- outfile:

  File path for saving/resuming the full edge-list trace (qs2 format).

- diagfile:

  Optional file path for saving convergence diagnostics (RDS format).

- diag:

  Logical; whether to compute diagnostics (currently unused, diagnostics
  are always computed).

- max_iter:

  Integer; total number of MCMC iterations per chain (ignored when
  `conv_thres` is set).

- nchains:

  Integer; number of independent chains.

- ncores:

  Integer; number of threads for C++ MCMC sampling.

- ncores_qs:

  Integer; number of threads for qs2 serialization.

- batch_size:

  Integer; number of iterations per batch.

- conv_thres:

  Numeric or `NULL`; if set, run until the ASDSF drops below this
  threshold instead of using `max_iter`.

- resume:

  Logical; if `TRUE`, resume from existing `outfile`.

## Value

A list of edge-list chains (one list of edge matrices per chain).
