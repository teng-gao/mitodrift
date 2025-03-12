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
#' @export 
optimize_tree_cpp = function(
    tree_init, logP, logA, max_iter = 100, ncores = 1, trace = TRUE, outfile = NULL, trace_interval = 5
) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    tree_init = reorder_phylo(tree_init)

    tree_current = tree_init
    max_current = sum(score_tree_bp_wrapper(tree_current$edge, logP, logA))
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

        scores = nni_cpp_parallel(tree_current$edge, logP, logA)
        
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

#' @export
reorder_phylo = function(phy) {
    phy_new = rlang::duplicate(phy, shallow = FALSE)
    phy_new$edge = reorderRcpp(phy$edge) %>% matrix(ncol = 2)
    return(phy_new)
}

#' @export
get_leaf_liks = function(mut_dat, vafs, ncores = 1) {

    variants = unique(mut_dat$variant)

    liks = mut_dat %>% 
        mutate(variant = factor(variant, variants)) %>%
        tidyr::complete(variant, cell, fill = list(vaf = 0, d = 0, a = 0)) %>%
        group_by(cell, variant) %>%
        group_modify(
            function(x, key) {
                l = sapply(vafs,
                    function(v) {
                        dbinom(x = x$a, size = x$d, prob = v)
                })
                tibble(l, vaf = vafs)
        }) %>%
        ungroup()

    liks = liks %>%
        split(.$variant) %>%
        mclapply(
            mc.cores = ncores,
            function(V) {
                V %>% reshape2::dcast(vaf ~ cell, value.var = 'l', fill = 1) %>%
                tibble::column_to_rownames('vaf') %>%
                as.matrix
            }
        )

    return(liks)
}

#' @export 
convert_liks_to_logP_list <- function(liks, phy) {
    
    E <- reorder_phylo(phy)$edge
    phy$node.label <- NULL
    
    P_all <- lapply(liks, function(liks_mut) {
        n_tips <- length(phy$tip.label)
        n_nodes <- phy$Nnode
        root_node <- E[nrow(E), 1]
        k <- nrow(liks_mut)
        
        P <- matrix(nrow = k, ncol = n_tips + n_nodes)
        rownames(P) <- rownames(liks_mut)

        # Tip likelihoods, internal node likelihoods, root node set up
        P[, 1:n_tips] <- liks_mut[, phy$tip.label]
        P[, (n_tips + 1):(n_tips + n_nodes)] <- 1/k
        P[, root_node] <- c(1, rep(0, k - 1))
        
        return(P)
    })

    logP_list <- lapply(P_all, function(P) {
        as.vector(log(t(P)))
    })
    
    return(logP_list)
}

#' @export 
run_tree_mcmc_cpp = function(phy, logP_list, A, max_iter = 100, nchains = 1, ncores = 1) {

    # edge_list = tree_mcmc_cpp(phy$edge, logP_list, t(log(A)), max_iter = max_iter)

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    edge_lists = tree_mcmc_parallel(phy$edge, logP_list, t(log(A)), max_iter = max_iter, nchains = nchains)

    res = edge_lists %>%
        lapply(function(elist) {
                lapply(
                    elist,
                    function(edges){
                        phy_new = rlang::duplicate(phy, shallow = FALSE)
                        E_new = matrix(edges, ncol = 2)
                        phy_new$edge = E_new
                    return(phy_new)
                })
            }
        )
    
    return(res)
}