# MitoDrift

# Installation
For now:
```R
devtools::install_local('.')
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

## One Step Wrapper
# run the complete MitoDrift pipeline
```
Rscript "$repo_dir/inst/bin/run_mitodrift_em.R" \
    --mut_dat "$mutfile" \
    --outdir "$OUTDIR" \
    --ncores "$NCORES" \
    --ncores_annot "$NCORES_ANNOT" \
    --k "$K" \
    --npop "$NPOP" \
    --eps "$EPS" \
    --err "$ERR" \
    --ngen "$NGEN" \
    --fit_params "$FIT_PARAMS" \
    --fit_param_max_iter "$FIT_PARAM_MAX_ITER" \
    --fit_param_epsilon "$FIT_PARAM_EPSILON" \
    --ml_iter "$ML_ITER" \
    --tree_mcmc_iter "$TREE_MCMC_ITER" \
    --tree_mcmc_chains "$TREE_MCMC_CHAINS" \
    --tree_mcmc_burnin "$TREE_MCMC_BURNIN" \
    --resume "$RESUME"
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

