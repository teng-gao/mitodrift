# Reorder a phylo object to postorder

Creates a deep copy of the phylogeny and reorders its edge matrix to
postorder using the compiled C++ helper `reorderRcpp`.

## Usage

``` r
reorder_phylo(phy)
```

## Arguments

- phy:

  A `phylo` object.

## Value

A new `phylo` object with edges in postorder.
