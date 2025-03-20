#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(igraph)
library(dplyr)
library(ggplot2)
library(ape)
library(tidygraph)
library(data.table)
library(stringr)
library(parallel)
library(tidytree)
library(glue)
library(patchwork)
library(phangorn)
library(ggtree)
library(ggraph)
library(mitodrift)
library(optparse)

# repo_dir = '/broad/sankaranlab/tgao/mitodrift/mitodrift'
# devtools::load_all(repo_dir)

option_list <- list(
    make_option(
        c("-m", "--mut_dat"),
        type = "character",
        default = NULL,
        help = "Mutation data file",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-i", "--max_iter"),
        type = "integer",
        default = 1000,
        help = "Maximum number of iteration",
        metavar = "INTEGER"
    ),
    make_option(
        c("-o", "--outfile"),
        type = "character",
        default = NULL,
        help = "Output file .RDS",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-p", "--ncores"),
        type = "integer",
        default = 1,
        help = "Number of cores to use",
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
    ),
    make_option(
        c("-f", "--freq_dat"),
        type = "character",
        default = NULL,
        help = "Mutation frequency table",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-r", "--resume"),
        type = "logical",
        default = FALSE,
        help = "Whether to resume from output file",
        metavar = "LOGICAL"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

init_method = 'nj'

for (arg in names(opts)) {
  message(paste0(arg, ": ", opts[[arg]]))
}

mut_dat = fread(opts$mut_dat) %>% select(variant, cell, d, a, vaf)

# Filter singletons
mut_dat = mut_dat %>% group_by(variant) %>% filter(sum(a>0)>1) %>% ungroup()

n_muts = length(unique(mut_dat$variant))
n_cells = length(unique(mut_dat$cell))

message(paste0('Number of mutations: ', n_muts))
message(paste0('Number of cells: ', n_cells))

if (n_cells == 0 | n_muts == 0) {
    stop('Zero cells or no shared mutations')
}

if (!opts$resume) {

    k = 20
    set.seed(0)
    A = mitodrift::get_transition_mat_wf(k = k, eps = opts$eps, N = opts$n_pop, n_gen = opts$n_gen)
    liks = mitodrift::get_leaf_liks(mut_dat, mitodrift::get_vaf_bins(k = k), eps = opts$seq_err)

    message('Initial clustering')

    # initial clustering
    vmat = mut_dat %>% 
        reshape2::dcast(variant ~ cell, fill = 0, value.var = 'vaf') %>% 
        tibble::column_to_rownames('variant')

    vmat[,'outgroup'] = 0

    dist_mat = vmat %>% as.matrix %>% t %>% dist(method = 'manhattan')

    if (init_method == 'hc') {
        
        phy_init = hclust(dist_mat, method = 'ward.D2') %>%
            as.phylo %>% root(outgroup = 'outgroup') %>% drop.tip('outgroup')

    } else if (init_method == 'nj') {

        phy_init = nj(dist_mat) %>%
            as.phylo %>% root(outgroup = 'outgroup') %>% drop.tip('outgroup')

    } else if (init_method == 'random') {

        phy_init = hclust(dist_mat, method = 'ward.D2') %>%
            as.phylo %>% root(outgroup = 'outgroup') %>% drop.tip('outgroup')

        set.seed(0)
        phyr = rtree(phy_init$Nnode+1)
        phyr$tip.label = sample(phy_init$tip.label)
        phyr$node.label = phy_init$node.label
        phy_init = phyr

    }

    message('Searching tree')
    if (is.null(opts$freq_dat) | opts$freq_dat == "") {
        logA_vec = as.vector(t(log(A)))
    } else {
        message('Variable mutation rate mode')
        mut_freq = fread(opts$freq_dat)

        eps_dict = mut_dat %>%
            distinct(variant) %>%
            left_join(mut_freq, by = join_by(variant)) %>%
            mutate(freq = ifelse(is.na(freq), 0, freq)) %>%
            mutate(eps_scaled = scale_eps(opts$eps, freq, mu = 2)) %>%
            {setNames(.$eps_scaled, .$variant)}

        logA_vec = lapply(
            names(liks),
            function(var) {
                AA = modify_A(A, eps = eps_dict[var])
                as.vector(t(log(AA)))
            }
        )
    }

    logP_list = convert_liks_to_logP_list(liks, phy_init)

    res = optimize_tree_cpp(
        phy_init, logP_list, logA_vec, ncores = opts$ncores,
        max_iter = opts$max_iter,
        outfile = opts$outfile,
        resume = opts$resume)

} else {

    res = optimize_tree_cpp(
        ncores = opts$ncores,
        max_iter = opts$max_iter,
        outfile = opts$outfile,
        resume = opts$resume)

}

saveRDS(res, opts$outfile)