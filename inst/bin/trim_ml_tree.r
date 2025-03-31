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
library(optparse)
library(mitodrift)

repo_dir = '/broad/sankaranlab/tgao/mitodrift/mitodrift'
R.utils::sourceDirectory(glue('{repo_dir}/R'))

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
        help = "Output file .RDS",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-b", "--burnin"),
        type = "integer",
        default = 0,
        help = "Output file .RDS",
        metavar = "INTEGER"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


for (arg in names(opts)) {
    message(paste0(arg, ": ", opts[[arg]]))
}

res_ml = readRDS(opts$tree_file)
tree_ml = res_ml$tree_list %>% .[[length(.)]]

res_mcmc = readRDS(opts$mcmc_file)
trees_mcmc = collect_chains(res_mcmc, burnin = opts$burnin)

tree_ml = add_clade_freq(tree_ml, trees_mcmc)

confs = seq(0, 1, by = 0.05)

trees_trim = lapply(
    confs,
    function(conf) {
        tree = trim_tree(tree_ml, conf = conf)
        tree$conf = conf
        return(tree)
    }
)

saveRDS(trees_trim, opts$outfile)