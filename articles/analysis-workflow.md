# Analysis workflow

This article walks through the standard post-inference analysis
workflow: **annotated tree → diagnostic PR curve → confidence trimming →
clone assignment → visualization**.

## Core concepts for interpretation

- **Initial tree topology**: a point-estimate starting tree constructed
  using neighbor joining (NJ) on continuous VAF matrices. This provides
  a fully-resolved (binary) initialization that empirically captures
  strong lineage signal before posterior sampling.
- **Posterior clade support**: per-node support values in
  `tree$node.label` (0–1) estimated from MCMC topology sampling.
- **Confidence-based topology refinement**: collapse internal edges
  below a support cutoff `τ` to obtain a refined lineage tree.

## Setup

``` r
library(mitodrift)
library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
```

We use a small in vitro LARRY barcode sample (200 cells, 186 variants)
bundled with the package.

``` r
data(pL1000_tree_annot)
data(pL1000_mut_dat)
```

In a typical workflow you would load the annotated tree from MitoDrift
output:

``` r
# tree_annot <- ape::read.tree("tree_annotated.newick")
```

## Visualize the full binary tree

[`plot_phylo_heatmap2()`](https://teng-gao.github.io/mitodrift/reference/plot_phylo_heatmap2.md)
displays the tree alongside a variant heteroplasmy heatmap. Setting
`node_conf = TRUE` colours internal nodes by their confidence score.

``` r
plot_phylo_heatmap2(
  pL1000_tree_annot,
  pL1000_mut_dat,
  node_conf = TRUE,
  branch_length = FALSE,
  title = "Full annotated tree"
)
```

![](analysis-workflow_files/figure-html/full-tree-1.png)

## Diagnostic: variant precision–recall curve

[`compute_variant_pr_curve()`](https://teng-gao.github.io/mitodrift/reference/compute_variant_pr_curve.md)
compares variant-defined cell partitions against tree clades across a
sweep of confidence cutoffs. This helps identify a threshold that
balances precision (are the clades real?) and recall (are we keeping
enough structure?).

``` r
pr_df <- compute_variant_pr_curve(pL1000_tree_annot, pL1000_mut_dat)
plot_prec_recall_vs_conf(pr_df, cutoff = 0.2)
```

![](analysis-workflow_files/figure-html/pr-curve-1.png)

## Trim tree

Collapse low-confidence nodes below the chosen threshold with
[`trim_tree()`](https://teng-gao.github.io/mitodrift/reference/trim_tree.md).

``` r
tree_trim <- trim_tree(pL1000_tree_annot, conf = 0.2)
```

Visualize the trimmed tree — polytomies replace poorly supported splits.

``` r
plot_phylo_heatmap2(
  tree_trim,
  pL1000_mut_dat,
  node_conf = TRUE,
  branch_length = FALSE,
  title = "Trimmed tree (conf >= 0.2)"
)
```

![](analysis-workflow_files/figure-html/trimmed-tree-1.png)

## Clone assignment

[`assign_clones_polytomy()`](https://teng-gao.github.io/mitodrift/reference/assign_clones_polytomy.md)
partitions tips into clones based on the polytomy structure of the
trimmed tree.

``` r
clone_df <- assign_clones_polytomy(tree_trim)
head(clone_df)
```

    ## # A tibble: 6 × 6
    ##   cell               clade clade_node annot  size  frac
    ##   <chr>              <chr>      <int> <chr> <int> <dbl>
    ## 1 CAACTAATCATTGACA-1 1            220 1        15 0.075
    ## 2 GTTCATTTCGGTTTGG-1 1            220 1        15 0.075
    ## 3 GGGCAATAGGCCCAGT-1 1            220 1        15 0.075
    ## 4 ATGTAAGCAATTGCGC-1 1            220 1        15 0.075
    ## 5 AGCTGCTCATAATTGC-1 1            220 1        15 0.075
    ## 6 GTTTCAGCAACACTTG-1 1            220 1        15 0.075

## Visualize clones

Colour cells by clone assignment on both a rectangular heatmap view and
a circular layout.

``` r
plot_phylo_heatmap2(
  tree_trim,
  pL1000_mut_dat,
  cell_annot = clone_df,
  branch_length = FALSE,
  title = "Clones on trimmed tree"
)
```

![](analysis-workflow_files/figure-html/clone-heatmap-1.png)

``` r
plot_phylo_circ(
  tree_trim,
  cell_annot = clone_df,
  title = "Circular layout"
)
```

![](analysis-workflow_files/figure-html/clone-circ-1.png)
