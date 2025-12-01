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
    library(mitodrift)
    library(optparse)
})

# repo_dir = '/broad/sankaranlab/tgao/mitodrift/mitodrift'
# R.utils::sourceDirectory(glue('{repo_dir}/R'))

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
        c("-n", "--npop"),
        type = "integer",
        default = 600,
        help = "Population size",
        metavar = "INTEGER"
    ),
    make_option(
        c("-g", "--ngen"),
        type = "integer",
        default = 100,
        help = "Number of generations",
        metavar = "INTEGER"
    ),
    make_option(
        c("-k", "--k"),
        type = "integer",
        default = 20,
        help = "Number of VAF bins",
        metavar = "INTEGER"
    ),
    make_option(
        c("-f", "--freq_dat"),
        type = "character",
        default = '',
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

mut_dat = fread(opts$mut_dat) %>% select(variant, cell, d, a) %>%
    mutate(vaf = a/d)

# Filter singletons
mut_dat = mut_dat %>% group_by(variant) %>% filter(sum(a>0)>1) %>% ungroup()

# Convert to matrix format
amat = long_to_mat(mut_dat, "a")
dmat = long_to_mat(mut_dat, "d")
vmat = long_to_mat(mut_dat, "vaf")

n_muts = nrow(amat)
n_cells = ncol(amat)

message(paste0('Number of mutations: ', n_muts))
message(paste0('Number of cells: ', n_cells))

if (n_cells == 0 | n_muts == 0) {
    stop('Zero cells or no shared mutations')
}

if (!opts$resume) {

    outdir <- dirname(opts$outfile)
    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    nj_ncores <- suppressWarnings(as.integer(opts$ncores))
    if (is.na(nj_ncores) || nj_ncores < 1L) {
        nj_ncores <- 1L
    }
    core_label <- if (nj_ncores > 1L) 'cores' else 'core'
    message(glue('Building initial tree using NJ ({nj_ncores} {core_label})'))
    phy_init = make_rooted_nj(vmat, ncores = nj_ncores)

    A = get_transition_mat_wf_hmm(k = opts$k, eps = opts$eps, N = opts$npop, ngen = opts$ngen)
    liks = get_leaf_liks_mat(amat, dmat, get_vaf_bins(k = opts$k), eps = opts$seq_err, log = TRUE)

    message('Searching for ML tree')
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

    logP_list = convert_logliks_to_logP_list(liks, phy_init)

    res = optimize_tree_cpp(
        phy_init, logP_list, logA_vec, 
        ncores = opts$ncores,
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

message('All done!')

saveRDS(res, opts$outfile)