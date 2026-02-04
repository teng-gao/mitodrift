# Plot Circular Phylogenetic Tree with Annotations

Creates a circular (fan) layout phylogenetic tree plot with optional
annotations, feature heatmaps, and node confidence values.

## Usage

``` r
plot_phylo_circ(
  gtree,
  cell_annot = NULL,
  tip_annot = NULL,
  feature_mat = NULL,
  branch_scores = NULL,
  node_conf = FALSE,
  conf_label = FALSE,
  title = "",
  pwidth_annot = 0.25,
  pwidth_feature = 0.25,
  branch_width = 0.3,
  dot_size = 1,
  conf_min = 0,
  conf_max = 0.5,
  annot_pal = NULL,
  offset = 0.05,
  width = 0.8,
  label_size = 2,
  rescale = FALSE,
  limits = c(-2, 2),
  feature_legend_title = "Score",
  flip = TRUE,
  annot_legend = TRUE,
  feature_legend = TRUE,
  layered = FALSE,
  smooth_k = 0,
  ladderize = TRUE,
  open_angle = 0,
  feature_axis_label = TRUE,
  feature_axis_angle = 30
)
```

## Arguments

- gtree:

  A `phylo` object or `tbl_graph` object representing the tree.

- cell_annot:

  A data frame or list of data frames containing cell annotations.
  Columns should include 'cell' and 'annot'.

- tip_annot:

  A data frame containing tip annotations. Columns should include 'cell'
  and 'annot'.

- feature_mat:

  A matrix of features to plot as a heatmap (rows are features, columns
  are cells).

- branch_scores:

  Optional branch scores (currently unused).

- node_conf:

  Logical. If `TRUE` and `gtree` is a `tbl_graph`, plots node
  confidence.

- conf_label:

  Logical. If `TRUE`, adds text labels for node confidence.

- title:

  Character string for the plot title.

- pwidth_annot:

  Numeric. Width of the annotation ring(s).

- pwidth_feature:

  Numeric. Width of the feature heatmap ring.

- branch_width:

  Numeric. Line width for tree branches.

- dot_size:

  Numeric. Size of tip/node points.

- conf_min:

  Numeric. Minimum value for confidence color scale.

- conf_max:

  Numeric. Maximum value for confidence color scale.

- annot_pal:

  A vector or list of vectors specifying colors for annotations.

- offset:

  Numeric. Offset distance between tree and annotation rings.

- width:

  Numeric. Width of heatmap cells.

- label_size:

  Numeric. Size of text labels.

- rescale:

  Logical. If `TRUE`, rescales features (z-score) for the heatmap.

- limits:

  Numeric vector of length 2. Limits for the feature heatmap color
  scale.

- feature_legend_title:

  Character string. Title for the feature heatmap legend.

- flip:

  Logical. If `TRUE`, flips the tree direction.

- layered:

  Logical. If `TRUE`, plots annotations as layered tiles.

- smooth_k:

  Integer. Number of neighbors for smoothing features. 0 means no
  smoothing.

- ladderize:

  Logical. If `TRUE`, ladderizes the tree.

- open_angle:

  Numeric. Angle of the opening in the fan layout.

- feature_axis_label:

  Logical. If `TRUE`, show feature heatmap axis labels; hide when
  `FALSE`.

- feature_axis_angle:

  Numeric. Rotation angle (degrees) for feature axis text.

- legend:

  Logical. If `TRUE`, shows legends.

## Value

A `ggtree` plot object.
