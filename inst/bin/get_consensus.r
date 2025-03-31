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
        c("-i", "--input"),
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
    ),
    make_option(
        c("-p", "--ncores"),
        type = "integer",
        default = 1,
        help = "Number of cores to use",
        metavar = "INTEGER"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


for (arg in names(opts)) {
  message(paste0(arg, ": ", opts[[arg]]))
}

res_mcmc = readRDS(opts$input)
trees_mcmc = collect_chains(res_mcmc, burnin = opts$burnin)

ps = c(0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)

gtrees = mclapply(
        ps,
        mc.cores = opts$ncores,
        function(p) {
            mitodrift::get_consensus(trees_mcmc, p = p)
    }) %>% setNames(ps)

saveRDS(gtrees, opts$outfile)