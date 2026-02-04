# Collapse weak clades below a confidence threshold.

Given a fully binary phylogeny with per-node confidence scores stored in
`tree$node.label`, this helper collapses every internal node whose
confidence falls below `conf`, effectively pruning low-support clades
while retaining higher-confidence structure.

## Usage

``` r
trim_tree(tree, conf, collapse_trivial = TRUE)
```

## Arguments

- tree:

  A binary `phylo` object whose `node.label` vector stores
  posterior/confidence values for internal nodes.

- conf:

  Numeric threshold in \$0,1\$. Internal nodes with confidence below
  this value are collapsed.

- collapse_trivial:

  Logical; if `TRUE`, collapses the trivial root-adjacent singleton
  split (1 tip vs the rest) into a root polytomy.

## Value

A renumbered `phylo` object with low-confidence nodes collapsed.
