# Package index

## Core workflow

- [`MitoDrift`](https://sankaranlab.github.io/mitodrift/reference/MitoDrift.md)
  : MitoDrift R6 Class

## Clades and support

- [`add_clade_freq()`](https://sankaranlab.github.io/mitodrift/reference/add_clade_freq.md)
  : Add clade frequencies to a phylogenetic tree
- [`collect_chains()`](https://sankaranlab.github.io/mitodrift/reference/collect_chains.md)
  : Collect MCMC chains into a multiPhylo object

## Topology refinement and clone assignment

- [`trim_tree()`](https://sankaranlab.github.io/mitodrift/reference/trim_tree.md)
  : Collapse weak clades below a confidence threshold.
- [`assign_clones_polytomy()`](https://sankaranlab.github.io/mitodrift/reference/assign_clones_polytomy.md)
  : Assign clone IDs to a tree allowing small polytomies.

## Diagnostics

- [`compute_variant_pr_curve()`](https://sankaranlab.github.io/mitodrift/reference/compute_variant_pr_curve.md)
  : Compute variant-to-tree precision/recall/F1 curve across confidence
  cutoffs
- [`plot_prec_recall_vs_conf()`](https://sankaranlab.github.io/mitodrift/reference/plot_prec_recall_vs_conf.md)
  : Plot precision/recall/F1 versus confidence cutoff

## Visualization

- [`plot_phylo_circ()`](https://sankaranlab.github.io/mitodrift/reference/plot_phylo_circ.md)
  : Plot Circular Phylogenetic Tree with Annotations
- [`plot_phylo_heatmap2()`](https://sankaranlab.github.io/mitodrift/reference/plot_phylo_heatmap2.md)
  : Plot a phylogenetic tree with VAF heatmap and annotations
- [`make_clade_pal()`](https://sankaranlab.github.io/mitodrift/reference/make_clade_pal.md)
  : Build a named Brewer palette with optional cycling behavior

## Data utilities

- [`mat_to_long()`](https://sankaranlab.github.io/mitodrift/reference/mat_to_long.md)
  : Convert allele count matrix (amat) and total count matrix (dmat) to
  long format

- [`long_to_mat()`](https://sankaranlab.github.io/mitodrift/reference/long_to_mat.md)
  : Convert a long format mutation data frame to a matrix

- [`phylo_to_gtree()`](https://sankaranlab.github.io/mitodrift/reference/phylo_to_gtree.md)
  :

  Convert a `phylo` object to a `tbl_graph`

- [`add_node_names()`](https://sankaranlab.github.io/mitodrift/reference/add_node_names.md)
  :

  Add "Node" labels to the internal nodes of a phylo tree

## Example data

- [`pL1000_tree_annot`](https://sankaranlab.github.io/mitodrift/reference/pL1000_tree_annot.md)
  : Annotated phylogenetic tree for pL1000_200_1
- [`pL1000_mut_dat`](https://sankaranlab.github.io/mitodrift/reference/pL1000_mut_dat.md)
  : Mutation count data for pL1000_200_1
