# Convert a `phylo` object to a `tbl_graph`

Builds a `tbl_graph` with tip nodes labeled by `phy$tip.label` and
internal nodes labeled from `phy$node.label`. If internal labels are
missing, they are created via
[`add_node_names()`](https://teng-gao.github.io/mitodrift/reference/add_node_names.md).

## Usage

``` r
phylo_to_gtree(phy)
```

## Arguments

- phy:

  A `phylo` object.

## Value

A `tbl_graph` with `nodes` and `edges` from the input tree.
