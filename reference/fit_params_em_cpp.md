# Fit tree parameters using EM (BP-backed, native ordering)

A faster EM that uses compute_node_edge_beliefs_bp2() directly. It
relies on the inherent ordering used by the C++ routine:

- node_beliefs\[locus\] is an n x C matrix in APE node-id order (1..n)

- edge_beliefs\[locus\] is an m x C x C array where the m edges
  correspond to the post-ordered edge list internal to the C++ routine
  (order is irrelevant because we sum over edges).

## Usage

``` r
fit_params_em_cpp(
  tree_fit,
  amat,
  dmat,
  initial_params = c(ngen = 100, log_eps = log(0.001), log_err = log(0.001)),
  lower_bounds = c(ngen = 1, log_eps = log(1e-12), log_err = log(1e-12)),
  upper_bounds = c(ngen = 1000, log_eps = log(0.2), log_err = log(0.2)),
  max_iter = 10,
  k = 20,
  npop = 600,
  ncores = 1,
  epsilon = 0.001,
  trace = TRUE
)
```

## Arguments

- tree_fit:

  phylogenetic tree (will be renumbered)

- amat:

  alternative allele counts

- dmat:

  total depth matrix

- initial_params:

  named numeric: ngen, log_eps, log_err

- lower_bounds:

  named numeric lower bounds (same names)

- upper_bounds:

  named numeric upper bounds (same names)

- max_iter:

  maximum EM iterations

- k:

  number of hidden states (VAF bins)

- npop:

  population size for WF-HMM transition

- ncores:

  kept for API compatibility (unused here)

- epsilon:

  convergence threshold on parameter deltas

- trace:

  if TRUE, returns trace dataframe

## Value

either named vector of final params (ngen, eps, err) or list(par=...,
trace=...)
