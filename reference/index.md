# Package index

## Core workflow

- [`MitoDrift`](https://teng-gao.github.io/mitodrift/reference/MitoDrift.md)
  : MitoDrift R6 Class
- [`fit_params_em_cpp()`](https://teng-gao.github.io/mitodrift/reference/fit_params_em_cpp.md)
  : Fit tree parameters using EM (BP-backed, native ordering)
- [`decode_tree()`](https://teng-gao.github.io/mitodrift/reference/decode_tree.md)
  : Decode a tree using CRF belief propagation (R version)
- [`make_rooted_nj()`](https://teng-gao.github.io/mitodrift/reference/make_rooted_nj.md)
  : Make a rooted NJ tree
- [`optimize_tree_cpp()`](https://teng-gao.github.io/mitodrift/reference/optimize_tree_cpp.md)
  : Optimize tree topology using C++ NNI moves
- [`run_tree_mcmc_batch()`](https://teng-gao.github.io/mitodrift/reference/run_tree_mcmc_batch.md)
  : Run tree-topology MCMC in batches with convergence monitoring

## Clades and support

- [`compute_target_tree_asdsf()`](https://teng-gao.github.io/mitodrift/reference/compute_target_tree_asdsf.md)
  : Compute ASDSF across chains using clades from a target tree

- [`add_clade_freq()`](https://teng-gao.github.io/mitodrift/reference/add_clade_freq.md)
  : Add clade frequencies to a phylogenetic tree

- [`add_node_names()`](https://teng-gao.github.io/mitodrift/reference/add_node_names.md)
  :

  Add "Node" labels to the internal nodes of a phylo tree

- [`collect_chains()`](https://teng-gao.github.io/mitodrift/reference/collect_chains.md)
  : Collect MCMC chains into a multiPhylo object

## Trimming and clones

- [`trim_tree()`](https://teng-gao.github.io/mitodrift/reference/trim_tree.md)
  : Collapse weak clades below a confidence threshold.
- [`trim_tree_exp()`](https://teng-gao.github.io/mitodrift/reference/trim_tree_exp.md)
  : Collapse branches with high expected mis-assignments.
- [`assign_clones_polytomy()`](https://teng-gao.github.io/mitodrift/reference/assign_clones_polytomy.md)
  : Assign clone IDs to a tree allowing small polytomies.

## Likelihood helpers

- [`get_leaf_liks_mat_cpp()`](https://teng-gao.github.io/mitodrift/reference/get_leaf_liks_mat_cpp.md)
  : Leaf likelihood matrix via C++ backend
- [`get_transition_mat_wf_hmm()`](https://teng-gao.github.io/mitodrift/reference/get_transition_mat_wf_hmm.md)
  : Get the transition matrix for WF model with HMM (with caching) TODO:
  add log option for small probabilities
- [`convert_logliks_to_logP_list()`](https://teng-gao.github.io/mitodrift/reference/convert_logliks_to_logP_list.md)
  : Convert log-likelihood matrices to a log-probability list
  (row-major)
- [`get_vaf_bins()`](https://teng-gao.github.io/mitodrift/reference/get_vaf_bins.md)
  : Generate VAF bin midpoints
- [`logSumExp()`](https://teng-gao.github.io/mitodrift/reference/logSumExp.md)
  : logSumExp function for a vector

## Visualization

- [`plot_phylo_circ()`](https://teng-gao.github.io/mitodrift/reference/plot_phylo_circ.md)
  : Plot Circular Phylogenetic Tree with Annotations
- [`plot_phylo_heatmap2()`](https://teng-gao.github.io/mitodrift/reference/plot_phylo_heatmap2.md)
  : Plot a phylogenetic tree with VAF heatmap and annotations

## Data utilities

- [`mat_to_long()`](https://teng-gao.github.io/mitodrift/reference/mat_to_long.md)
  : Convert allele count matrix (amat) and total count matrix (dmat) to
  long format

- [`long_to_mat()`](https://teng-gao.github.io/mitodrift/reference/long_to_mat.md)
  : Convert a long format mutation data frame to a matrix

- [`phylo_to_gtree()`](https://teng-gao.github.io/mitodrift/reference/phylo_to_gtree.md)
  :

  Convert a `phylo` object to a `tbl_graph`

- [`reorder_phylo()`](https://teng-gao.github.io/mitodrift/reference/reorder_phylo.md)
  : Reorder a phylo object to postorder
