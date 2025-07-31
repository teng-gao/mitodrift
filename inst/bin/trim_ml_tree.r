#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages({
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
    library(optparse)
    library(mitodrift)
})
# repo_dir = '/broad/sankaranlab/tgao/mitodrift/mitodrift'
# R.utils::sourceDirectory(glue('{repo_dir}/R'))

option_list <- list(
    make_option(
        c("-t", "--tree_file"),
        type = "character",
        help = "ML result RDS file",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-c", "--mcmc_file"),
        type = "character",
        help = "MCMC result RDS file",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-o", "--outfile"),
        type = "character",
        help = "Output tree file .newick",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-b", "--burnin"),
        type = "integer",
        default = 0,
        help = "Burnin",
        metavar = "INTEGER"
    ),
    make_option(
        c("-m", "--max_iter"),
        type = "integer",
        default = 1e8,
        help = "Maximum iteration",
        metavar = "INTEGER"
    ),
    make_option(
        c("-u", "--use_nj"),
        type = "logical",
        default = FALSE,
        help = "Whether to use the NJ tree as the initial tree",
        metavar = "LOGICAL"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

for (arg in names(opts)) {
    message(paste0(arg, ": ", opts[[arg]]))
}

res_ml = readRDS(opts$tree_file)

if (opts$use_nj) {
    tree = res_ml$tree_list[[1]]
} else {
    tree = res_ml$tree_list %>% .[[length(.)]]
}

res_mcmc = readRDS(opts$mcmc_file)
trees_mcmc = collect_chains(res_mcmc, burnin = opts$burnin, max_iter = opts$max_iter)
tree = add_clade_freq(tree, trees_mcmc)

write.tree(tree, opts$outfile)