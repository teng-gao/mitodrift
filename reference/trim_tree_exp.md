# Collapse branches with high expected mis-assignments.

Computes the expected number of incorrectly assigned tips for each
internal branch (based on (1 - p) \* clade size, where p is the branch
confidence) and collapses branches whose expectation exceeds a tolerance
derived from the total number of tips.

## Usage

``` r
trim_tree_exp(tree, tol)
```

## Arguments

- tree:

  A binary `phylo` object with internal node confidence scores in
  `tree$node.label`.

- tol:

  Numeric tolerance expressed as a fraction of total tips; the threshold
  is `length(tree$tip.label) * tol` expected errors.

## Value

A renumbered `phylo` object with overly uncertain branches collapsed.
