# Assign clone IDs to a tree allowing small polytomies.

Assign clone IDs to a tree allowing small polytomies.

## Usage

``` r
assign_clones_polytomy(tree, k = Inf, paraphyletic = FALSE, return_df = TRUE)
```

## Arguments

- tree:

  Rooted `phylo` object with bifurcating or polytomous structure.

- k:

  Positive scalar; maximum clade size that collapses into one clone.

- paraphyletic:

  Logical; if `TRUE`, leftover unassigned tips below a node are grouped
  into one clone, otherwise each becomes a singleton.

- return_df:

  Logical; if `TRUE`, returns a data frame annotated per tip, otherwise
  an integer vector of clone assignments.

## Value

Either a data frame with columns `cell`, `clade`, `annot`, `size`,
`frac` or an integer vector of clone IDs indexed by tip order.

## Details

Finds the root (node never used as a child), recurses over internal
nodes, and assigns clone IDs whenever a clade has `<= k` tips. Tips
directly attached to the root are always singleton clones. Ensures all
tips receive an assignment, warning if late singletons are required.
