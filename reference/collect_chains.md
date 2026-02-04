# Collect MCMC chains into a multiPhylo object

Reconstructs `phylo` trees from raw edge-list chains, applies burn-in
removal and iteration truncation, then pools all chains into a single
`multiPhylo` object.

## Usage

``` r
collect_chains(edge_list_all, phy_init, burnin = 0, max_iter = Inf)
```

## Arguments

- edge_list_all:

  List of chains, each a list of 2-column integer edge matrices.

- phy_init:

  The initial `phylo` object whose tip labels and metadata are used to
  reconstruct full `phylo` objects.

- burnin:

  Integer; number of initial samples to discard from each chain.

- max_iter:

  Numeric; maximum iteration to retain (samples beyond this are
  dropped).

## Value

A `multiPhylo` object containing the pooled post-burn-in trees.
