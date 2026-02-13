# Plot precision/recall/F1 versus confidence cutoff

Build a line plot of precision, recall, and F1 across confidence cutoffs
from a PR-curve data frame (for example, output from
[`compute_variant_pr_curve()`](https://teng-gao.github.io/mitodrift/reference/compute_variant_pr_curve.md),
`compute_tree_pr_curve()`, or `compute_pr_from_jaccard_mat()`).

## Usage

``` r
plot_prec_recall_vs_conf(
  pr_df,
  sample_name = NULL,
  cutoff = NULL,
  x_axis_scale = c("log10", "linear"),
  log10_min = 1e-06,
  j_thres = NULL,
  legend = TRUE
)
```

## Arguments

- pr_df:

  Data frame containing at least columns `conf_cutoff`, `precision`,
  `recall`, and `f1`. If present, `fraction_well_matched` is used for
  optional annotation.

- sample_name:

  Optional title string shown above the panel.

- cutoff:

  Optional numeric cutoff to highlight with a vertical dashed line and
  label near the top of the panel.

- x_axis_scale:

  X-axis scale: `"log10"` (default) or `"linear"`.

- log10_min:

  Lower bound used for log10 plotting so zero cutoffs can be displayed
  safely.

- j_thres:

  Optional Jaccard threshold used only for display annotation with
  `fraction_well_matched`.

- legend:

  Logical; if `TRUE`, show metric legend on the right.

## Value

A `ggplot2` object.
