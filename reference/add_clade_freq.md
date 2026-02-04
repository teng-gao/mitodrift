# Add clade frequencies to a phylogenetic tree

Computes the frequency (support) of each clade in a reference phylogeny
(`phy`) across a list of phylogenetic trees (`edge_list`), and adds
these frequencies as node labels to the reference tree.

## Usage

``` r
add_clade_freq(phy, edge_list, rooted = TRUE, ncores = 1)
```

## Arguments

- phy:

  A reference phylogeny of class `phylo`.

- edge_list:

  A list of edge matrices (or phylogenetic trees) to compare against the
  reference tree.

- rooted:

  Logical; whether to treat the trees as rooted. Default is `TRUE`.

- ncores:

  Integer; number of cores to use for parallel computation. Default is
  `1` (no parallelization).

## Value

A phylo object with clade frequencies added as node labels.
