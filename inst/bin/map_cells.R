#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(stringr)
    library(dplyr)
    library(parallel)
    library(ape)
    library(glue)
    library(tidygraph)
    library(optparse)
    library(TreeTools)
    library(mitodrift)
})

repo_dir = "/lab-share/Hem-Sankaran-e2/Public/projects/tgao/tools/mitodrift"
R.utils::sourceDirectory(glue('{repo_dir}/R'))

option_list <- list(
    make_option(
        c("-b", "--mitodrift_obj"),
        type = "character",
        help = "Mitodrift object",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-a", "--amat"),
        type = "character",
        help = "Allele count matrix",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-d", "--dmat"),
        type = "character",
        help = "Total count matrix",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-x", "--prefix"),
        type = "character",
        help = "Prefix for cell names in the mitodrift object",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-m", "--mut_dat"),
        type = "character",
        help = "Mutation data file (alternative to providing amat/dmat separately)",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-c", "--cell_list"),
        type = "character",
        help = "List of cells to map; by default all cells in the count matrices will be mapped",
        default = NULL,
        metavar = "CHARACTER"
    ),
    make_option(
        c("-o", "--outfile"),
        type = "character",
        help = "Output file .CSV (cell, tip, prob)",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-p", "--ncores"),
        type = "integer",
        help = "Number of cores",
        metavar = "INTEGER"
    ),
    make_option(
        c("-l", "--leaf_only"),
        type = "logical",
        help = "Whether to only map to leaf nodes",
        default = TRUE,
        metavar = "LOGICAL"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

for (arg in names(opts)) {
    if (opts[[arg]] == "" | opts[[arg]] == 'NULL') {
        opts[[arg]] = NULL
    }
    message(paste0(arg, ": ", opts[[arg]]))
}

outfile = opts$outfile
ncores = opts$ncores
leaf_only = opts$leaf_only

if (!is.null(opts$mut_dat)) {
    mut_dat = fread(opts$mut_dat)
    amat = long_to_mat(mut_dat, 'a')
    dmat = long_to_mat(mut_dat, 'd')
} else if (!is.null(opts$amat) & !is.null(opts$dmat)) {
    amat = fread(opts$amat) %>% tibble::column_to_rownames('Variants') %>% as.matrix()
    dmat = fread(opts$dmat) %>% tibble::column_to_rownames('Variants') %>% as.matrix()
} else {
    stop('Either mut_dat or amat/dmat must be provided')
}

md = readRDS(opts$mitodrift_obj)

if (!is.null(opts$prefix)) {
    md$tree_annot$tip.label = paste0(opts$prefix, '_', md$tree_annot$tip.label)
    colnames(md$amat) = paste0(opts$prefix, '_', colnames(md$amat))
    colnames(md$dmat) = paste0(opts$prefix, '_', colnames(md$dmat))
}

if (!all(colnames(md$amat) %in% colnames(amat)) | !all(rownames(md$amat) %in% rownames(amat))) {
    stop('Allele count matrices does not match mitodrift object')
}

if (!is.null(opts$cell_list)) {
    cells_map_all = readLines(opts$cell_list)
    message(glue('read {length(cells_map_all)} cells from cell list'))
} else {
    cells_map_all = colnames(amat)
}

tree_guide = md$tree_annot
cells_guide = tree_guide$tip.label
cells_map_all = cells_map_all[!cells_map_all %in% cells_guide]
cells_map_all = cells_map_all %>% intersect(colnames(amat))

# fetch model parameters
logA_vec = md$logA
k = md$model_params[['k']]
err = md$model_params[['err']]
vaf_bins = get_vaf_bins(k = k)

message('filtering cells by shared mutations')

# filter cells with shared variants
variants_guide = rownames(md$amat)

variants_new = amat[,cells_map_all] %>% apply(2, function(x){names(x[x>0])})

n_shared = sapply(
    variants_new,
    function(vars) {
        variants_share = intersect(variants_guide, vars)
        length(variants_share)
    }
)

cells_map = names(n_shared[n_shared >= 1])
n_cells = length(cells_map)

amat = amat[,c(cells_guide, cells_map)] %>% .[rowSums(.) > 1,,drop=FALSE]
dmat = dmat[rownames(amat),colnames(amat)]

message(glue('Conducting mapping of {n_cells} cells onto a guide tree with {length(cells_guide)} cells'))
message(glue('Number of variants in guide tree: {length(variants_guide)}'))
message(glue('Number of variants in full matrix: {nrow(amat)}'))

logliks_all = get_leaf_liks_mat_cpp(amat, dmat, vaf_bins, log = TRUE, eps = err)

message(glue('Starting mapping'))

res_all = lapply(
        cells_map,
        function(cell_map) {
    
            message(cell_map)
            
            # subset likelihood by relevant cells and variants
            cells_comb = c(tree_guide$tip.label, cell_map)
            variants_comb = union(variants_new[[cell_map]], variants_guide) %>% intersect(names(logliks_all))
            liks_comb = logliks_all[variants_comb] %>% lapply(function(liks){liks[,colnames(liks) %in% cells_comb,drop=FALSE]})
            
            probs = map_cell_to_tree(
                tree_guide,
                cell_map,
                liks_comb, 
                logA_vec,
                ncores = ncores,
                leaf_only = leaf_only
            )
    
            res = tibble(p = unname(probs), guide_node = names(probs), cell_map = cell_map)

            fwrite(res, outfile, append = TRUE)

            return(res)
    
        }
    ) %>%
    bind_rows()

message('Done!')
fwrite(res_all, outfile)

