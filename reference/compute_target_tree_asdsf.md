# Compute ASDSF across chains using clades from a target tree

Compute ASDSF across chains using clades from a target tree

## Usage

``` r
compute_target_tree_asdsf(
  phy_target,
  edge_list_chains,
  rooted = TRUE,
  normalize = TRUE,
  ncores = 1,
  min_freq = 0
)
```

## Arguments

- phy_target:

  A rooted `phylo` object defining the reference clades. Assumed to be
  in postorder.

- edge_list_chains:

  List of chains, each a list of edge matrices (2-column integer
  matrices).

- rooted:

  Logical; treat trees as rooted when matching clades. Default `TRUE`.

- normalize:

  Logical; pass through to `prop_clades_par` (default `TRUE`).

- min_freq:

  Minimum clade frequency threshold used when averaging SDs. Default
  `0.1`.

## Value

A list with elements `asdsf`, `per_clade_sd`, `keep_mask`, and
`freq_matrix`.
