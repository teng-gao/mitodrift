# MitoDrift

MitoDrift reconstructs single-cell lineage relationships from mitochondrial DNA (mtDNA) mutations by modeling heteroplasmy drift and measurement noise {using a WF-HMM model}. It produces **confidence-calibrated phylogenies** and downstream summaries (trimmed trees, clade support, and clone partitions). The input to MitoDrift are mtDNA allele counts from single-cell genomics data {add list of protocols, e.g. mtscATAC-seq, ReDeeM, etc}, and outputs are a phylogenetic tree with posterior clade support, a confidence-trimmed tree, and clone-level summaries/visualizations {this is repetitive}

---

## Installation

```r
# from a local checkout
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_local(".")
```

For faster MCMC trace IO, install **qs2** with TBB:

```r
remotes::install_cran("qs2", type = "source", configure.args = "--with-TBB --with-simd=AVX2")
```

---

## Quick start

This uses a tiny example dataset packaged under `inst/extdata/`.

{why are these in comments? should be shell codeblocks}
# You can run this from a terminal in the package root:
#
#   Rscript inst/bin/run_mitodrift_em.R \
#     --mut_dat inst/extdata/small_test_mut_dat.csv \
#     --outdir mitodrift_demo \
#     --tree_mcmc_iter 50 \
#     --tree_mcmc_chains 2 \
#     --tree_mcmc_burnin 10
#
# Outputs include:
#   - mitodrift_object.rds
#   - tree files and diagnostics

Once the run finishes:

```r
library(mitodrift)
mut_dat <- read.csv(
  system.file("extdata", "small_test_mut_dat.csv", package = "mitodrift")
)
md <- readRDS("mitodrift_demo/mitodrift_object.rds")
phy_trim <- trim_tree(md$tree_annot, conf = 0.5)

plot_phylo_heatmap2(
  phy_trim,
  mut_dat,
  dot_size = 1,
  branch_width = 0.3,
  branch_length = FALSE,
  node_conf = TRUE,
  het_max = 1
)
```


## Input data requirements

### Long format (`mut_dat`)
One row per cell–variant pair, with counts for alternate allele (`a`) and total depth (`d`).

| cell      | variant       | a    | d    |
|-----------|---------------|------|------|
| PD45534aj | MT_10448_T_C  | 0    | 2348 |
| PD45534aj | MT_11787_T_C  | 1462 | 2348 |

**Important**: Include rows where `a = 0` so that **every** cell × variant combination is represented (the observation model uses total depth).

### Matrix format (`amat`, `dmat`)
Wide matrices with `variant` in the first column and cells in subsequent columns. Use the same layout for `amat` (alt counts) and `dmat` (total depth).

Helper conversion:

```r
long_dat <- mat_to_long(amat, dmat)
```


## Inference settings & diagnostics

Key MCMC settings (from `run_mitodrift_em.R`):
- `tree_mcmc_burnin`: number of initial samples to discard.
- `tree_mcmc_chains`: number of independent runs (10 typical; 50 for thorough analysis).
- `tree_mcmc_batch_size`: how often to check convergence.

Convergence diagnostics:
- ASDSF (`conv_thres`) summarizes agreement across chains.
- Diagnostics are written to `tree_mcmc_diag.rds`.
- For interrupted runs, use `--resume TRUE` to continue from saved traces.


---

## Core concepts for interpretation
- **Initial tree topology**: {Constructed using NJ on continuous VAF matrices, we find that this approach recovers the most lienage signal for a full binary tree. this is the point estimate for tree toplogy}
- **Posterior clade support**: per-node support values in `tree$node.label` (0–1) estimated from MCMC topology sampling.
- **Confidence trimming**: collapse internal edges below a support cutoff `τ` to obtain a confidence-calibrated lineage tree.

---

---

## Clone assignment workflow

1) **Trim** low-confidence edges:
```r
phy_trim <- trim_tree(phy, conf = τ)
```

2) **Optional**: subset or drop root singletons (if needed for a specific analysis).

3) **Assign clones** (top-level root-descending clades):
```r
clade_df <- assign_clones_polytomy(phy_trim, k = Inf, return_df = TRUE)
```

**Interpretation**: root singletons reflect unresolved relationships and should not be interpreted as biologically unique clones.

---

---

## References

Check out our manuscript (preprint/published version will be added here) 
