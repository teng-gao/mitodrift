# Compute variant-to-tree precision/recall/F1 curve across confidence cutoffs

Construct a variant VAF matrix from long-format mutation counts, compare
variant partitions against tree clades using Jaccard similarity, and
sweep node-confidence cutoffs to summarize precision, recall, and F1.

## Usage

``` r
compute_variant_pr_curve(
  phy_annot,
  mut_dat,
  min_vaf = 0.01,
  min_cells = 2,
  j_thres = 0.66,
  n_points = 50,
  ncores = 1
)
```

## Arguments

- phy_annot:

  A rooted `phylo` object with clade confidence stored in `node.label`.

- mut_dat:

  Long-format mutation table containing at least columns `cell`,
  `variant`, `a`, and `d`.

- min_vaf:

  Minimum VAF used when defining variant-positive cells (default
  `0.01`).

- min_cells:

  Minimum number of variant-positive cells required for a variant to be
  evaluated (default `2`).

- j_thres:

  Jaccard threshold for declaring a match between a variant split and a
  tree clade (default `0.66`).

- n_points:

  Number of confidence cutoffs to evaluate (default `50`).

- ncores:

  Number of cores for sweeping confidence cutoffs. Uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) when
  `ncores > 1` on non-Windows systems.

## Value

A data frame with columns `conf_cutoff`, `n_pred`, `n_tp`, `n_tru`,
`n_recover`, `precision`, `recall`, `f1`, and `fraction_well_matched`.
