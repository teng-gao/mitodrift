# MitoDrift

# Installation
For now:
```R
devtools::install_local('.')
```
Make sure you install [qs2](https://github.com/qsbase/qs2) with TBB binding to enable faster MCMC trace file saving/reading:
```
remotes::install_cran("qs2", type = "source", configure.args = "--with-TBB --with-simd=AVX2")
```

# Usage
## Input format
1. `mut_dat` includes the below columns, where `a` is the alternate allele count and `d` is the total depth.

| cell      | variant       | a    | d    |
|-----------|---------------|------|------|
| PD45534aj | MT_10448_T_C  | 0    | 2348 |
| PD45534aj | MT_11787_T_C  | 1462 | 2348 |
| PD45534aj | MT_11790_T_C  | 0    | 2338 |
| PD45534aj | MT_12242_A_G  | 0    | 2277 |
| PD45534aj | MT_12630_G_A  | 0    | 2382 |
| PD45534aj | MT_15791_A_G  | 0    | 2459 |

Note that depth information should be provided even if there is no alternate allele detected at the variant position (a=0), so that there is a row for each cell x variant combination.

2. `amat` / `dmat` inputs are wide matrices (CSV). The first column must be the variant identifier, subsequent columns represent cells. Each row corresponds to one variant:

| variant     | cell_1 | cell_2 | cell_3 |
|-------------|--------|--------|--------|
| MT_10448_T_C| 0      | 0      | 3      |
| MT_11787_T_C| 12     | 8      | 15     |
| MT_11790_T_C| 0      | 1      | 0      |

Use the same layout for `dmat`, replacing counts with total depths.

# One Step Wrapper


```
inst/bin/run_mitodrift_em.R
Usage: run_mitodrift_em.R [options]


Options:
	-m CHARACTER, --mut_dat=CHARACTER
		Mutation table (CSV with columns: variant, cell, d, a)

	-A CHARACTER, --amat=CHARACTER
		Alternative allele count matrix (CSV; first column = variant)

	-D CHARACTER, --dmat=CHARACTER
		Total depth matrix (CSV; first column = variant)

	-o CHARACTER, --outdir=CHARACTER
		Output directory for results (default: mitodrift_results)

	-p INTEGER, --ncores=INTEGER
		Number of cores for ML/MCMC (default: 1)

	-q INTEGER, --ncores_annot=INTEGER
		Number of cores for branch confidence annotation (default: 1)

	-u INTEGER, --ncores_em=INTEGER
		Number of cores for EM parameter fitting (default: 1)

	-k INTEGER, --k=INTEGER
		Number of VAF bins (default: 20)

	-n INTEGER, --npop=INTEGER
		Population size (default: 600)

	-e DOUBLE, --eps=DOUBLE
		Mutation rate per branch (default: 0.001)

	-s DOUBLE, --err=DOUBLE
		Sequencing error rate (default: 0)

	-g INTEGER, --ngen=INTEGER
		Number of generations (default: 100)

	-f LOGICAL, --fit_params=LOGICAL
		Fit parameters using EM (default: TRUE)

	-t INTEGER, --fit_param_max_iter=INTEGER
		Maximum EM iterations (default: 10)

	-c DOUBLE, --fit_param_epsilon=DOUBLE
		EM convergence threshold (default: 1e-3)

	-i INTEGER, --ml_iter=INTEGER
		Maximum iterations for tree optimization (default: 100)

	-j INTEGER, --tree_mcmc_iter=INTEGER
		Maximum iterations for phylogenetic MCMC (default: 100)

	-l INTEGER, --tree_mcmc_chains=INTEGER
		Number of MCMC chains (default: 1)

	-b INTEGER, --tree_mcmc_burnin=INTEGER
		Burn-in iterations for phylogenetic MCMC (default: 0)

	-d INTEGER, --tree_mcmc_batch_size=INTEGER
		Batch size for phylogenetic MCMC (default: 1000)

	-y LOGICAL, --tree_mcmc_diag=LOGICAL
		Run diagnostics (e.g., ASDSF) after each MCMC batch (default: TRUE)

	-r LOGICAL, --resume=LOGICAL
		Resume from existing outputs (default: FALSE)

	-h, --help
		Show this help message and exit

```

```
Rscript "$repo_dir/inst/bin/run_mitodrift_em.R" \
    --mut_dat "$mutfile" \  # or replace with: --amat "$amat" --dmat "$dmat"
    --outdir "$outdir" \
    --ncores "$ncores" \
    --ncores_em "$ncores_em" \
    --ncores_annot "$ncores_annot" \
    --k "$k" \
    --npop "$npop" \
    --eps "$eps" \
    --err "$err" \
    --ngen "$ngen" \
    --fit_params "$fit_params" \
    --fit_param_max_iter "$fit_param_max_iter" \
    --fit_param_epsilon "$fit_param_epsilon" \
    --ml_iter "$ml_iter" \
    --tree_mcmc_iter "$tree_mcmc_iter" \
    --tree_mcmc_chains "$tree_mcmc_chains" \
    --tree_mcmc_burnin "$tree_mcmc_burnin" \
    --tree_mcmc_batch_size "$tree_mcmc_batch_size" \
    --tree_mcmc_diag "$tree_mcmc_diag" \
    --resume "$resume"
```

# Visualizing results
## ML tree
```
# KX003_1_tree_list is an example output of ML tree estimation module
tree_list = KX003_1_tree_list$tree_list
tree_ml = tree_list %>% .[[length(.)]]

plot_phylo_heatmap2(
    tree_ml, 
    KX003_1_mut_dat,
    dot_size = 0.1,
    branch_width = 0.3,
    branch_length = F,
    het_max = 1,
    title = 'ML'
)
```
Result:
![image](https://github.com/user-attachments/assets/c06daeeb-2dd4-4f35-bb8a-6e14c295f80a)

## Trimmed (consensus) ML tree with clade confidence
```
# Saved MitoDrift object in the output folder
md = readRDS('{outdir}/mitodrift_object.rds')

phy_trim = trim_tree(md$tree_annot, conf = 0.5)

plot_phylo_heatmap2(
    phy_trim, # choose tree with branch confidence cutoff = 0.5
    mut_dat, # input mutation dataframe for the ML module
    dot_size = 1,
    branch_width = 0.3,
    branch_length = F,
    node_conf = T, 
    het_max = 1,
    title = 'Trimmed ML'
)
```
Result:
![image](https://github.com/user-attachments/assets/62834948-2227-41f1-9826-b0080bad7b0d)

