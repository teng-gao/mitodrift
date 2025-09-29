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
    library(fmcmc)
    library(CRF)
    library(optimParallel)
})

# repo_dir = "/lab-share/Hem-Sankaran-e2/Public/projects/tgao/tools/mitodrift"
# devtools::load_all(repo_dir)

# Define command line options
option_list <- list(
    make_option(
        c("-m", "--mut_dat"),
        type = "character",
        help = "Mutation data file (CSV format with columns: variant, cell, d, a)",
        metavar = "CHARACTER"
    ),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        default = "mitodrift_results",
        help = "Output directory for results",
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
        c("-q", "--ncores_annot"),
        type = "integer",
        default = 1,
        help = "Number of cores to use for branch confidence annotation",
        metavar = "INTEGER"
    ),
    make_option(
        c("-u", "--ncores_em"),
        type = "integer",
        default = 1,
        help = "Number of cores to use for EM parameter fitting",
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
        c("-n", "--npop"),
        type = "integer",
        default = 600,
        help = "Population size",
        metavar = "INTEGER"
    ),
    make_option(
        c("-e", "--eps"),
        type = "double",
        default = 0.001,
        help = "Mutation rate per branch",
        metavar = "DOUBLE"
    ),
    make_option(
        c("-s", "--err"),
        type = "double",
        default = 0,
        help = "Sequencing error rate",
        metavar = "DOUBLE"
    ),
    make_option(
        c("-g", "--ngen"),
        type = "integer",
        default = 100,
        help = "Number of generations",
        metavar = "INTEGER"
    ),
    make_option(
        c("-f", "--fit_params"),
        type = "logical",
        default = TRUE,
        help = "Whether to fit parameters using EM",
        metavar = "LOGICAL"
    ),
    make_option(
        c("-t", "--fit_param_max_iter"),
        type = "integer",
        default = 10,
        help = "Maximum EM iterations for parameter fitting",
        metavar = "INTEGER"
    ),
    make_option(
        c("-c", "--fit_param_epsilon"),
        type = "double",
        default = 1e-3,
        help = "Convergence threshold for EM parameter fitting",
        metavar = "DOUBLE"
    ),
    make_option(
        c("-i", "--ml_iter"),
        type = "integer",
        default = 100,
        help = "Maximum iterations for tree optimization",
        metavar = "INTEGER"
    ),
    make_option(
        c("-j", "--tree_mcmc_iter"),
        type = "integer",
        default = 100,
        help = "Maximum iterations for phylogenetic MCMC",
        metavar = "INTEGER"
    ),
    make_option(
        c("-l", "--tree_mcmc_chains"),
        type = "integer",
        default = 1,
        help = "Number of MCMC chains for phylogenetic sampling",
        metavar = "INTEGER"
    ),
    make_option(
        c("-b", "--tree_mcmc_burnin"),
        type = "integer",
        default = 0,
        help = "Burnin for phylogenetic MCMC",
        metavar = "INTEGER"
    ),
    make_option(
        c("-d", "--tree_mcmc_batch_size"),
        type = "integer",
        default = 1000,
        help = "Batch size for phylogenetic MCMC",
        metavar = "INTEGER"
    ),
    make_option(
        c("-y", "--tree_mcmc_diag"),
        type = "logical",
        default = TRUE,
        help = "Whether to run diagnostics (e.g., ASDSF) after each MCMC batch",
        metavar = "LOGICAL"
    ),
    make_option(
        c("-r", "--resume"),
        type = "logical",
        default = FALSE,
        help = "Whether to resume from existing files",
        metavar = "LOGICAL"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Print parameters
message("=== MitoDrift Analysis Parameters ===")
for (arg in names(opts)) {
    message(paste0(arg, ": ", opts[[arg]]))
}
if (opts$resume) {
    message("RESUME MODE: Will attempt to resume from existing files")
}

# Check required parameters
if (is.null(opts$mut_dat)) {
    stop("Mutation data file must be provided with -m/--mut_dat")
}

# Create output directory
if (!dir.exists(opts$outdir)) {
    dir.create(opts$outdir, recursive = TRUE)
    message("Created output directory: ", opts$outdir)
}

# Define output file paths
ml_trace_file <- file.path(opts$outdir, "tree_ml_trace.rds")
mcmc_trace_file <- file.path(opts$outdir, "tree_mcmc_trace.rds")
annot_tree_file <- file.path(opts$outdir, "tree_annotated.newick")
mitodrift_object_file <- file.path(opts$outdir, "mitodrift_object.rds")

# Step 1: Load mutation data
message("\n=== Creating MitoDrift object ===")
mut_dat <- fread(opts$mut_dat) %>% 
    group_by(variant) %>%
    filter(sum(a > 0) > 1) %>%
    ungroup()
message("Number of unique variants: ", length(unique(mut_dat$variant)))
message("Number of unique cells: ", length(unique(mut_dat$cell)))

# Step 2: Create MitoDrift object
if (opts$resume) {
    message("=== Resuming from existing MitoDrift object ===")
    md <- readRDS(mitodrift_object_file)
    
    # Check if model components are initialized
    if (is.null(md$logP) || is.null(md$logA)) {
        stop("Model components not initialized")
    }
} else {
    md <- MitoDrift$new(
        mut_dat = mut_dat,
        model_params = c(
            eps = opts$eps,
            err = opts$err,
            npop = opts$npop,
            ngen = opts$ngen,
            k = opts$k
        )
    )
    saveRDS(md, mitodrift_object_file)

    if (opts$fit_params) {
        message("\n=== Fitting parameters using EM ===")
        fitted_params <- md$fit_params_em(
            tree_fit = md$tree_init,
            max_iter = opts$fit_param_max_iter,
            epsilon = opts$fit_param_epsilon,
            ncores = opts$ncores_em
        )
    } else {
        message("\n=== Skipping parameter fitting ===")
    }

    message("\n=== Initializing model components ===")
    md$make_model(ncores = opts$ncores)
    saveRDS(md, mitodrift_object_file)
}

message("\n=== Optimizing tree ===")
md$optimize_tree(
    max_iter = opts$ml_iter,
    ncores = opts$ncores,
    outfile = ml_trace_file,
    resume = opts$resume
)

saveRDS(md, mitodrift_object_file)

message("\n=== Running phylogenetic MCMC ===")
md$run_mcmc(
    max_iter = opts$tree_mcmc_iter,
    nchains = opts$tree_mcmc_chains,
    ncores = opts$ncores,
    batch_size = opts$tree_mcmc_batch_size,
    diag = opts$tree_mcmc_diag,
    outfile = mcmc_trace_file,
    resume = opts$resume
)

saveRDS(md, mitodrift_object_file)

message("\n=== Annotating tree with clade frequencies ===")
md$annotate_tree(
    burnin = opts$tree_mcmc_burnin,
    ncores = opts$ncores_annot
)

write.tree(md$tree_annot, annot_tree_file)
saveRDS(md, mitodrift_object_file)

message("\n=== MitoDrift analysis completed successfully! ===")

