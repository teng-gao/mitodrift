#' @import dplyr
#' @import tidygraph
#' @import stringr
#' @importFrom igraph vcount ecount E V V<- E<- 
#' @importFrom phangorn upgma 
#' @importFrom ape root drop.tip nj
#' @importFrom parallelDist parDist
#' @importFrom stats na.omit reorder setNames
#' @useDynLib mitodrift
NULL

# to fix: apparently the init tree has to be rooted otherwise to_phylo_reoder won't work. 
optimize_tree_cpp = function(
    tree_init, logP, logA, max_iter = 100, ncores = 1, trace = TRUE, outfile = NULL, trace_interval = 5
) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    tree_init$edge = reorderRcpp(tree_init$edge)

    tree_current = tree_init
    max_current = -Inf
    tree_list = list()
    runtime = c(0,0,0)

    for (i in 1:max_iter) {

        ptm = proc.time()

        if (trace) {
            tree_list = c(tree_list, list(tree_current))
            if (!is.null(outfile)) {
                if (i == 1 | i %% trace_interval == 0) {
                    saveRDS(tree_list, outfile)
                }
            }
        }

        message(paste(i, round(max_current, 4), paste0('(', signif(unname(runtime[3]),2), 's', ')')))

        scores = nni_cpp_parallel(tree_current, logP, logA)
        
        if (max(scores) > max_current) {
            max_id = which.max(scores)
            if (max_id %% 2 == 0) {pair_id = 2} else {pair_id = 1}
            tree_current$edge = matrix(nnin_cpp(tree_current$edge, ceiling(max_id/2))[[pair_id]], ncol = 2)
            tree_current$logZ = max_current = max(scores)
        } else {
            break()
        }

        runtime = proc.time() - ptm
        
    }

    if (trace) {
        class(tree_list) = 'multiPhylo'
        return(tree_list)
    } else {
        return(tree_current)
    }
}

# to fix: apparently the init tree has to be rooted otherwise to_phylo_reoder won't work. 
optimize_tree_rcpp = function(
    tree_init, logP, logA, max_iter = 100, ncores = 1, trace = TRUE, outfile = NULL, trace_interval = 5
) {

    tree_init$edge = reorderRcpp(tree_init$edge)

    tree_current = tree_init
    max_current = -Inf
    tree_list = list()
    runtime = c(0,0,0)

    for (i in 1:max_iter) {

        ptm = proc.time()

        if (trace) {
            tree_list = c(tree_list, list(tree_current))
            if (!is.null(outfile)) {
                if (i == 1 | i %% trace_interval == 0) {
                    saveRDS(tree_list, outfile)
                }
            }
        }

        message(paste(i, round(max_current, 4), paste0('(', signif(unname(runtime[3]),2), 's', ')')))

        # scores = nni_cpp_parallel(tree_current, logP, logA)

        trees_nei = TreeSearch::NNI(tree_current, edgeToBreak = -1)

        scores = mclapply(
            trees_nei,
            mc.cores = ncores,
            function(tree) {
                score_tree_bp_wrapper(reorderRcpp(tree$edge), logP, logA) %>% sum
            }
        ) %>% unlist()
        
        if (max(scores) > max_current) {
            tree_current$edge = reorderRcpp(trees_nei[[which.max(scores)]]$edge)
            tree_current$logZ = max_current = max(scores)
        } else {
            break()
        }

        runtime = proc.time() - ptm
        
    }

    if (trace) {
        class(tree_list) = 'multiPhylo'
        return(tree_list)
    } else {
        return(tree_current)
    }
}
