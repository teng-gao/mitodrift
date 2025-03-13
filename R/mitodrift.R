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
get_leaf_liks = function(mut_dat, vafs, eps = 0, ncores = 1) {

    variants = unique(mut_dat$variant)

    liks = mut_dat %>% 
        mutate(variant = factor(variant, variants)) %>%
        tidyr::complete(variant, cell, fill = list(vaf = 0, d = 0, a = 0)) %>%
        split(.$variant) %>%
        mclapply(
            mc.cores = ncores,
            function(mut_dat_var) {
                mut_dat_var %>%
                group_by(cell, variant) %>%
                group_modify(
                    function(x, key) {
                        l = sapply(vafs,
                            function(v) {
                                dbinom(x = x$a, size = x$d, prob = pmin(v + eps, 1))
                        })
                        tibble(l, vaf = vafs)
                }) %>%
                ungroup()
        }) %>%
        bind_rows()
        
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
get_vaf_bins = function(k) {
    bins = seq(0, 1, 1/k)
    bins = c(0,bins,1)
    vafs = sapply(1:(length(bins)-1), function(i){(bins[i] + bins[i+1])/2})
    return(vafs)
}

#' @export 
get_transition_mat_wf = function(k, eps = 0.01, N = 100, n_rep = 1e4, n_gen = 100) {

    A = matrix(NA, ncol = k + 2, k + 2)
    bins = seq(0, 1, 1/k)
    bins = c(0,bins,1)

    for(i in 1:(length(bins)-1)) {

        p = mean(c(bins[i], bins[i+1]))

        p_gen = p
        for (gen in 1:n_gen) {
            x = rbinom(n_rep, N, p_gen)
            p_gen = x/N   
        }
        
        for(j in 1:(length(bins)-1)) {
            
            xstart = as.integer(N*bins[j]) + 1
            xend = as.integer(N*bins[j+1])
            xstart = min(xend, xstart)

            if (j == (length(bins)-2)) {
                xend = xend - 1
            }

            A[i, j] = length(x[x >= xstart & x <= xend])/n_rep

        }
    }
    
    A[1,] = c(1-eps, rep(eps/(ncol(A)-1), ncol(A)-1))
    A[nrow(A),] = rev(A[1,])

    if (n_gen == 0) {
        A = diag(k+2)
    }

    vs = sapply(1:(length(bins)-1), function(i){(bins[i] + bins[i+1])/2})
    colnames(A) = vs
    rownames(A) = vs

    return(A)

}

#' @export 
run_mitodrift = function(
    mut_dat, n_gen = 100, eps = 0.001, n_pop = 600, k = 20, seq_err = 0, max_iter = 100,
    init_method = 'nj', ncores = 1, outfile = NULL
) {

    set.seed(0)
    A = get_transition_mat_wf(k = k, eps = eps, N = n_pop, n_gen = n_gen)
    liks = get_leaf_liks(mut_dat, get_vaf_bins(k = k), eps = seq_err, ncores = ncores)

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

        phy_init = ape::nj(dist_mat) %>%
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
    logA_vec = as.vector(t(log(A)))
    logP_list = convert_liks_to_logP_list(liks, phy_init)

    tree_list = optimize_tree_cpp(
        phy_init, logP_list, logA_vec, ncores = ncores,
        max_iter = max_iter, 
        outfile = outfile)

    return(tree_list)
}

################################### MCMC ######################################


#' @export 
run_tree_mcmc_cpp = function(phy, logP_list, logA_vec, max_iter = 100, nchains = 1, ncores = 1) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)
    
    edge_lists = tree_mcmc_parallel(phy$edge, logP_list, logA_vec, max_iter = max_iter, nchains = nchains)

    res = edge_lists %>%
        lapply(function(elist) {
                lapply(
                    elist,
                    function(edges){
                    phy_new = attach_edges(phy, edges)
                    return(phy_new)
                })
            }
        )
    
    return(res)
}

attach_edges = function(phy, edges) {

    phy_new = rlang::duplicate(phy, shallow = FALSE)
    E_new = matrix(edges, ncol = 2)
    phy_new$edge = E_new

    return(phy_new)
}

# to fix: apparently the init tree has to be rooted otherwise to_phylo_reoder won't work. 
#' @export 
run_tree_mcmc_r = function(
    gtree_init, A, liks, max_iter = 100, nchains = 1, ncores = 1, trace = TRUE, outfile = NULL,
    move_type = 'NNI'
) {

    # Rearrangements have to be unrooted
    if (move_type == 'NNI') {
        propose_tree = TreeSearch::NNI
    } else {
        stop('Invalid move type')
    }

    res = mclapply(
        1:nchains,
        mc.cores = ncores,
        function(s) {

            set.seed(s)

            # need to change this to pre-allocation
            tree_list = vector("list", max_iter + 1)
            phy_list = vector("list", max_iter + 1)
            tree_list[[1]] = gtree_init
            phy_list[[1]] = as.phylo(gtree_init)

            for (i in 1:max_iter) {
                
                phy_new = propose_tree(phy_list[[i]])
                gtree_new = decode_tree(phy_new, A, liks, score_only = TRUE)
                
                l_0 = sum(tree_list[[i]]$logZ)
                l_1 = sum(gtree_new$logZ)
                probab = exp(l_1 - l_0)

                if (runif(1) < probab) {
                    tree_list[[i+1]] = gtree_new
                    phy_list[[i+1]] = phy_new
                } else {
                    tree_list[[i+1]] = tree_list[[i]]
                    phy_list[[i+1]] = phy_list[[i]]
                }
            }

            if (trace &!is.null(outfile)) {
                outdir = dirname(outfile)
                fname = basename(outfile)
                saveRDS(phy_list, glue('{outdir}/chain{s}_{fname}'))
            }

            return(list('tree_list' = tree_list, 'phy_list' = phy_list))
        }
    )

    return(res)
}

# to fix: apparently the init tree has to be rooted otherwise to_phylo_reoder won't work. 
#' @export
run_tree_mcmc = function(
    phy_init, logP_list, logA_vec, max_iter = 100, nchains = 1, ncores = 1, trace = TRUE, outfile = NULL
) {

    res = mclapply(
        1:nchains,
        mc.cores = ncores,
        function(s) {

            elist = tree_mcmc_cpp(phy_init$edge, logP_list, logA_vec, max_iter = max_iter, seed = s)

            tree_list = lapply(
                elist,
                function(edges){
                    attach_edges(phy_init, edges)
            })

            if (trace &!is.null(outfile)) {
                outdir = dirname(outfile)
                fname = basename(outfile)
                saveRDS(tree_list, glue('{outdir}/chain{s}_{fname}'))
            }

            return(tree_list)
        }
    )

    return(res)
}

#' @export
collect_chains = function(res, burnin = 0) {

    mcmc_trees = res %>% lapply(function(trees){
            trees = trees[(burnin+1):length(trees)]
            lapply(trees, function(tree){
                tree$edge.length = NULL
                tree$node.label = NULL
                tree$nodes = NULL
                tree$edge = TreeTools::RenumberTree(tree$edge[,1], tree$edge[,2])
                return(tree)
            })
        }) %>% unlist(recursive = F)

        return(mcmc_trees)
    
}

