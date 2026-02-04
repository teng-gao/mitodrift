# Plot a phylogenetic tree with VAF heatmap and annotations

Plot a phylogenetic tree with VAF heatmap and annotations

## Usage

``` r
plot_phylo_heatmap2(
  phylo,
  df_var = NULL,
  branch_width = 0.25,
  root_edge = TRUE,
  dot_size = 1,
  ylim = NULL,
  min_cells = 1,
  tip_annot = NULL,
  annot_scale = NULL,
  feature_mat = NULL,
  feature_limits = c(-2, 2),
  feature_scale = NULL,
  rescale = FALSE,
  title = NULL,
  ytitle = NULL,
  xtitle = NULL,
  label_site = FALSE,
  cell_annot = NULL,
  tip_lab = FALSE,
  node_lab = FALSE,
  layered = FALSE,
  annot_bar_height = 0.1,
  clade_bar_height = 1,
  feature_height = 1,
  het_max = 0.1,
  conf_min = 0,
  conf_max = 1,
  conf_label = FALSE,
  branch_length = TRUE,
  node_conf = FALSE,
  annot_pal = NULL,
  annot_legend = FALSE,
  label_group = FALSE,
  text_size = 3,
  annot_title_size = text_size,
  node_label_size = 1,
  mut = NULL,
  mark_low_cov = FALSE,
  facet_by_group = FALSE,
  flip = TRUE,
  ladderize = TRUE,
  node_scores = NULL,
  node_score_limits = c(-4, 4),
  variants_highlight = NULL,
  show_variant_names = TRUE,
  show_tree_y_axis = FALSE,
  feature_legend = TRUE,
  raster = FALSE,
  raster_dpi = 300
)
```

## Arguments

- df_var:

  Optional data frame of variant calls with columns `cell`, `variant`,
  and either `vaf` or `a`/`d` for computing VAF.

- branch_width:

  Numeric branch line width for the tree.

- root_edge:

  Logical; included for compatibility with earlier calls.

- dot_size:

  Numeric size for tip/node points.

- ylim:

  Optional y-axis limits (unused in current implementation).

- min_cells:

  Integer minimum number of cells with VAF \> 0 required to keep a
  variant.

- tip_annot:

  Optional data frame of tip annotations with columns `cell` and
  `annot`.

- annot_scale:

  Optional ggplot scale for annotation colors.

- feature_mat:

  Optional matrix of features (rows = features, columns = cells) to plot
  as a heatmap.

- feature_limits:

  Numeric vector of length 2 giving limits for feature heatmap color
  scale.

- feature_scale:

  Optional ggplot scale for feature heatmap colors.

- rescale:

  Logical; if `TRUE`, z-score features per row before plotting.

- title:

  Optional plot title.

- ytitle:

  Optional VAF heatmap y-axis title.

- xtitle:

  Optional VAF heatmap x-axis title.

- label_site:

  Logical; included for compatibility with earlier calls.

- cell_annot:

  Optional data frame or list of data frames for annotation bars.

- tip_lab:

  Logical; if `TRUE`, show tip labels on the tree.

- node_lab:

  Logical; if `TRUE`, show internal node labels on the tree.

- layered:

  Logical; if `TRUE`, render annotation bars as layered tiles.

- annot_bar_height:

  Numeric height for each annotation bar panel.

- clade_bar_height:

  Numeric height for clade bar panel (unused in current implementation).

- feature_height:

  Numeric height for the feature heatmap panel.

- het_max:

  Numeric maximum VAF for heatmap color scaling.

- conf_min:

  Numeric minimum value for confidence color scale.

- conf_max:

  Numeric maximum value for confidence color scale.

- conf_label:

  Logical; if `TRUE`, label node confidence values.

- branch_length:

  Logical; if `FALSE`, drop branch lengths before plotting.

- node_conf:

  Logical; if `TRUE`, plot node confidence (requires `tbl_graph`).

- annot_pal:

  Optional palette (vector or list) for annotation bars.

- annot_legend:

  Logical; if `TRUE`, show annotation legends.

- label_group:

  Logical; if `TRUE`, label groups in annotation bars.

- text_size:

  Numeric base text size for labels.

- annot_title_size:

  Numeric size for annotation strip titles.

- node_label_size:

  Numeric size for node/tip labels.

- mut:

  Optional mutation column name in `tree` to color nodes by VAF.

- mark_low_cov:

  Logical; if `TRUE`, mark low-coverage cells on the VAF heatmap.

- facet_by_group:

  Logical; if `TRUE`, facet VAF heatmap by `group` column.

- flip:

  Logical; if `TRUE`, flip the tree orientation.

- ladderize:

  Logical; if `TRUE`, ladderize the tree.

- node_scores:

  Optional named numeric vector of node/tip scores used to color tree
  branches. Names may be tip labels, internal node labels, or numeric
  node IDs.

- node_score_limits:

  Numeric vector of length 2 giving limits for node score colors.

- variants_highlight:

  Optional vector of variant names to bold on the heatmap axis.

- show_variant_names:

  Logical; if `TRUE`, show variant tick labels on the VAF heatmap; hide
  when `FALSE`.

- show_tree_y_axis:

  Logical; if `TRUE`, show y-axis ticks/labels on the tree panel.

- feature_legend:

  Logical; if `TRUE`, display the legend for the feature heatmap;
  otherwise hide it.

- raster:

  Logical; if `TRUE`, rasterize each plot panel via `ggrastr`.

- raster_dpi:

  Numeric DPI to use when rasterizing panels.

- tree:

  A `phylo` object or `tbl_graph` tree.

## Value

A patchwork/ggplot object combining the tree, annotations, and heatmaps.
