# Add "Node" labels to the internal nodes of a phylo tree

Add "Node" labels to the internal nodes of a phylo tree

## Usage

``` r
add_node_names(tree, prefix = "Node", start_from_tip = TRUE)
```

## Arguments

- tree:

  A phylo object

- prefix:

  Character prefix for node names (default "Node")

- start_from_tip:

  Logical; if `TRUE` (default), node numbering starts from `ntip + 1`
  (standard ape convention); if `FALSE`, numbering starts from 1.

## Value

The same phylo object, with tree\$node.label set to prefix + node
numbers
