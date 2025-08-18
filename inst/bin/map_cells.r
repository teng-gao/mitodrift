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
options(repr.matrix.max.cols=50, repr.matrix.max.rows=50)
devtools::load_all('/broad/sankaranlab/tgao/mitodrift/mitodrift')

home = "/broad/sankaranlab/tgao"

option_list <- list(
    make_option(
        c("-t", "--tree"),
        type = "character",
        help = "Tree in Newick format",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-m", "--mut_dat"),
        type = "character",
        help = "Mutation data file",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-l", "--cell_list"),
        type = "character",
        help = "List of cells to map",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-o", "--outfile"),
        type = "character",
        help = "Output file .RDS",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-p", "--ncores"),
        type = "integer",
        help = "Number of cores",
        metavar = "INTEGER"
    ),
    make_option(
        c("-e", "--eps"),
        type = "double",
        default = 1e-3,
        help = "mutation rate",
        metavar = "DOUBLE"
    ),
    make_option(
        c("-s", "--seq_err"),
        type = "double",
        default = 0,
        help = "sequencing error rate",
        metavar = "DOUBLE"
    ),
    make_option(
        c("-n", "--n_pop"),
        type = "integer",
        default = 600,
        help = "Population size",
        metavar = "INTEGER"
    ),
    make_option(
        c("-g", "--n_gen"),
        type = "integer",
        default = 100,
        help = "Number of generations",
        metavar = "INTEGER"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

for (arg in names(opts)) {
  message(paste0(arg, ": ", opts[[arg]]))
}

tree_guide = ape::read.tree(opts$tree)

k = 20
set.seed(0)
A = mitodrift::get_transition_mat_wf(k = k, eps = opts$eps, N = opts$n_pop, n_gen = opts$n_gen)
logA_vec = as.vector(t(log(A)))

cells_map_all = readLines(opts$cell_list)
outfile = opts$outfile
ncores = opts$ncores
mut_dat_all = fread(opts$mut_dat)

map_cell_to_tree = function(tree, cell_new, liks_map, logA_vec, ncores = 1) {
    
    # generate all possible assignment to tip
    trees_new = lapply(
        1:length(tree$tip.label),
        function(i) {
            TreeTools::AddTip(tree, where = i, label = cell_new)
        }
    )
    
    class(trees_new) = 'multiPhylo'

    # subset likelihoods
    logP_list_map = convert_logliks_to_logP_list(liks_map, trees_new[[1]])

    # score assignments
    scores = mclapply(
            trees_new,
            mc.cores = ncores,
            function(phy) {
                score_tree_bp_wrapper(reorderRcpp(phy$edge), logP_list_map, logA_vec)
            }
        ) %>% unlist()
    
    probs = exp(scores - logSumExp(scores)) %>% setNames(tree$tip.label)

    return(probs)
}

set.seed(0)

message('filtering cells by shared mutations')
cells_guide = tree_guide$tip.label

# do not map cells that are already in the guide tree
cells_map_all = cells_map_all[!cells_map_all %in% cells_guide]

# filter cells with shared variants
variants_guide = mut_dat_all %>%
    filter(cell %in% cells_guide) %>%
    filter(a>0) %>%
    pull(variant)

variants_new = mut_dat_all %>%
    filter(cell %in% cells_map_all) %>%
    filter(a>0) %>% split(.$cell) %>% 
    lapply(function(x){x$variant})

n_shared = sapply(
    variants_new,
    function(vars) {
        variants_share = intersect(variants_guide, vars)
        length(variants_share)
    }
)

cells_map = names(n_shared[n_shared >= 1])
n_cells = length(cells_map)

mut_dat_map = mut_dat_all %>% 
    filter(cell %in% c(cells_guide, cells_map)) %>% 
    group_by(variant) %>%
    filter(sum(a>0)>1) %>%
    ungroup()

n_variants = length(unique(mut_dat_map$variant))

tree_guide = tree_guide %>% keep.tip(cells_guide[cells_guide %in% mut_dat_map$cell])
n_guide = length(tree_guide$tip.label)

message(glue('Conducting mapping of {n_cells} cells onto a guide tree with {n_guide} cells'))

message(glue('Number of shared variants: {n_variants}'))

liks_all = get_leaf_liks(mut_dat_map, mitodrift::get_vaf_bins(k = 20), log = TRUE, ncores = ncores, eps = opts$seq_err)

message(glue('Starting mapping'))

res_all = mclapply(
        cells_map,
        mc.cores = ncores,
        function(cell_new) {
    
            message(cell_new)
            
            cells_comb = c(tree_guide$tip.label, cell_new)

            variants_comb = mut_dat_map %>% 
                filter(cell %in% c(tree_guide$tip.label, cell_new)) %>% 
                group_by(variant) %>%
                summarise(n_cells = sum(a>0)) %>%
                filter(n_cells > 1) %>%
                pull(variant)

            liks_comb = liks_all[variants_comb] %>% lapply(function(liks){liks[,colnames(liks) %in% cells_comb,drop=F]})
            
            probs = map_cell_to_tree(
                tree_guide,
                cell_new, 
                liks_comb, 
                logA_vec
            )
    
            res = tibble(p = unname(probs), tip = names(probs), cell_new = cell_new)

            fwrite(res, outfile, append = TRUE)

            return(res)
    
        }
    ) %>%
    bind_rows()

message('Done!')
fwrite(res_all, outfile)

