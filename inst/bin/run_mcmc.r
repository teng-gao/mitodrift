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

repo_dir = '/broad/sankaranlab/tgao/mitodrift/mitodrift'
devtools::load_all(repo_dir)

option_list <- list(
    make_option(
        c("-t", "--tree_file"),
        type = "character",
        help = "Tree list file from mitodrift ML module",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-o", "--outfile"),
        type = "character",
        help = "Output file .RDS",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-i", "--max_iter"),
        type = "integer",
        default = 10000,
        help = "Maximum number of iteration",
        metavar = "INTEGER"
    ),
    make_option(
        c("-c", "--nchains"),
        type = "integer",
        default = 1000,
        help = "Maximum number of iteration",
        metavar = "INTEGER"
    ),
    make_option(
        c("-p", "--ncores"),
        type = "integer",
        default = 1,
        help = "Number of cores to use",
        metavar = "INTEGER"
    ),
    make_option(
        c("-u", "--use_nj"),
        type = "logical",
        default = FALSE,
        help = "Whether to use the NJ tree as the initial tree",
        metavar = "LOGICAL"
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

for (arg in names(opts)) {
  message(paste0(arg, ": ", opts[[arg]]))
}

res_ml = readRDS(opts$tree_file)

if (opts$use_nj) {
    message('Using NJ tree as initial tree')
    index = 1
} else {
    message('Using ML tree as initial tree')
    index = length(res_ml$tree_list)
}

phy_init = res_ml$tree_list[[index]]

logA_vec = res_ml$params$logA
logP_list = res_ml$params$logP

message('Running MCMC')

res_mcmc = mitodrift::run_tree_mcmc(
    phy_init, logP_list, logA_vec,
    max_iter = opts$max_iter, nchains = opts$nchains,
    ncores = opts$ncores, outfile = opts$outfile, resume = opts$resume)

message('Done!')