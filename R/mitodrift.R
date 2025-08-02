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
    tree_init = NULL, logP = NULL, logA = NULL, max_iter = 100, 
    outfile = NULL, resume = FALSE,
    ncores = 1, trace_interval = 5
) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    if (resume) {
        if (is.null(outfile)) {
            stop('Outfile must be provided if resume = FALSE')
        }
        message('Resuming from saved outfile')
        res = readRDS(outfile)
        tree_current = res$tree_list %>% .[[length(.)]]
        max_current = tree_current$logZ
        start_iter = length(res$tree_list)
        logP = res$params$logP
        logA = res$params$logA
    } else {
        start_iter = 1
    }

    if (is.list(logA)) {
        score_tree_func = score_tree_bp_wrapper_multi
        nni_func = nni_cpp_parallel_multi
    } else {
        score_tree_func = score_tree_bp_wrapper
        nni_func = nni_cpp_parallel
    }

    if (!resume) {
        if (is.null(tree_init) || is.null(logP) || is.null(logA)) {
            stop("tree_init, logP, and logA must be provided when resume is FALSE")
        }
        tree_init = reorder_phylo(tree_init)
        tree_init$edge.length = NULL
        tree_current = tree_init
        max_current = sum(score_tree_func(tree_current$edge, logP, logA))
        tree_current$logZ = max_current
        res = list(
            'params' = list('logP' = logP, 'logA' = logA),
            'tree_list' = phytools::as.multiPhylo(tree_current)
        )
    }

    runtime = c(0,0,0)

    if (start_iter > max_iter) {
        message(paste0("Already completed ", start_iter - 1, " iterations. No new iterations to run for max_iter = ", max_iter, "."))
        return(res)
    }

    for (i in start_iter:max_iter) {

        ptm = proc.time()

        message(paste(i, round(max_current, 4), paste0('(', signif(unname(runtime[3]),2), 's', ')')))

        scores = nni_func(tree_current$edge, logP, logA)
        
        if (max(scores) > max_current) {
            max_id = which.max(scores)
            if (max_id %% 2 == 0) {pair_id = 2} else {pair_id = 1}
            tree_current$edge = matrix(nnin_cpp(tree_current$edge, ceiling(max_id/2))[[pair_id]], ncol = 2)
            tree_current$logZ = max_current = max(scores)
        } else {
            break()
        }

        res$tree_list = res$tree_list %>% c(phytools::as.multiPhylo(tree_current))

        if (!is.null(outfile)) {
            if (i == 1 | i %% trace_interval == 0) {
                saveRDS(res, outfile)
            }
        }

        runtime = proc.time() - ptm
        
    }

    saveRDS(res, outfile)

    return(res)
}

# R version of the function
# to fix: apparently the init tree has to be rooted otherwise to_phylo_reoder won't work. 
#' @export
optimize_tree = function(
    gtree_init, A, liks, max_iter = 100, ncores = 1, trace = TRUE, outfile = NULL,
    trace_interval = 5, polytomy = FALSE, d = Inf, method = 'nni_multi', resume = FALSE
) {

    if (resume) {
        if (is.null(outfile)) {
            stop('Outfile must be provided if resume = FALSE')
        }
        message('Resuming from saved outfile')
        tree_list = readRDS(outfile)
        gtree = tree_list %>% .[[length(.)]]
        score = sum(gtree$logZ)
        # liks = res$params$A
        # A = res$params$liks
    } else {
        gtree = gtree_init
        tree_list = list()
    }

    runtime = c(0,0,0)

    for (i in 1:max_iter) {

        ptm = proc.time()

        if (trace) {
            tree_list = c(tree_list, list(gtree))
            if (!is.null(outfile)) {
                if (i == 1 | i %% trace_interval == 0) {
                    saveRDS(tree_list, outfile)
                }
            }
        }

        score = sum(gtree$logZ)

        message(paste(i, round(score, 4), paste0('(', signif(unname(runtime[3]),2), 's', ')')))

        if (polytomy) {
            if (method == 'nni_multi') {
                nei = nni_multi(as.phylo(gtree), d = d)
            } else {
                nei = nnt(as.phylo(gtree))
            }
        } else {
            nei = TreeSearch::NNI(as.phylo(gtree), edgeToBreak = -1)
        }
        
        gtrees_nei = mclapply(
            nei,
            mc.cores = ncores,
            function(tn) {
                decode_tree(tn, A, liks, score_only = TRUE)
            }
        )

        scores_nei = sapply(
            gtrees_nei,
            function(gtree) {
                sum(gtree$logZ)
            }
        )

        if (max(scores_nei) > score) {
            gtree = gtrees_nei[[which.max(scores_nei)]]
        } else {
            break()
        }

        runtime = proc.time() - ptm
        
    }

    if (trace) {
        return(tree_list)
    } else {
        return(gtree)
    }
}


#' @export
reorder_phylo = function(phy) {
    phy_new = rlang::duplicate(phy, shallow = FALSE)
    phy_new$edge = reorderRcpp(phy$edge) %>% matrix(ncol = 2)
    return(phy_new)
}

#' @export
get_leaf_liks = function(mut_dat, vafs, eps = 0, ncores = 1, log = FALSE) {

    variants = unique(mut_dat$variant)

    liks = mut_dat %>% 
        mutate(vaf = a/d) %>%
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
                                dbinom(x = x$a, size = x$d, prob = pmin(v + eps, 1), log = log)
                        })
                        tibble(l, vaf = vafs)
                }) %>%
                ungroup()
        }) %>%
        bind_rows()

    if (log) {
        na_fill = 0
    } else {
        na_fill = 1
    }
        
    liks = liks %>%
        split(.$variant) %>%
        mclapply(
            mc.cores = ncores,
            function(V) {
                V %>% reshape2::dcast(vaf ~ cell, value.var = 'l', fill = na_fill) %>%
                tibble::column_to_rownames('vaf') %>%
                as.matrix
            }
        )

    return(liks)
}

#' get leaft likelihoods given mutation data in matrix format
#' @param amat matrix of allele counts
#' @param dmat matrix of total counts
#' @param vafs vector of VAFs
#' @param eps error rate
#' @param log whether to return log likelihoods
#' @return list of likelihood matrices
#' @export
get_leaf_liks_mat = function(amat, dmat, vafs, eps = 0, ncores = 1, log = FALSE) {

    variants = rownames(amat)

    liks = mclapply(
        variants,
        mc.cores = ncores,
        function(v) {
            m = sapply(
                    vafs,
                    function(vaf) {
                        dbinom(x = amat[v,], 
                            size = dmat[v,],
                            prob = pmin(vaf + eps, 1),
                            log = log)
                    }
                )
            colnames(m) = vafs
            return(t(m))
        }) %>% setNames(variants)
    
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
        as.vector(t(log(P)))
    })
    
    return(logP_list)
}

#' @export 
convert_logliks_to_logP_list <- function(logliks, phy) {
    
    E <- reorder_phylo(phy)$edge
    phy$node.label <- NULL
    
    P_all <- lapply(logliks, function(liks_mut) {
        n_tips <- length(phy$tip.label)
        n_nodes <- phy$Nnode
        root_node <- E[nrow(E), 1]
        k <- nrow(liks_mut)
        
        P <- matrix(nrow = k, ncol = n_tips + n_nodes)
        rownames(P) <- rownames(liks_mut)

        # Tip likelihoods, internal node likelihoods, root node set up
        P[, 1:n_tips] <- liks_mut[, phy$tip.label]
        P[, (n_tips + 1):(n_tips + n_nodes)] <- log(1/k)
        P[, root_node] <- log(c(1, rep(0, k - 1)))
        
        return(P)
    })

    logP_list <- lapply(P_all, function(P) {
        as.vector(t(P))
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

#' Get the transition matrix for WF model with HMM
#' @param k number of VAF bins
#' @param eps error rate
#' @param N population size
#' @param n_gen number of generations
#' @return transition matrix
#' @export
get_transition_mat_wf_hmm <- function(k, eps, N, n_gen) {
    # For zero generations, the transition matrix is the identity matrix
    if (n_gen == 0) {
        A <- diag(k + 2)
    } else {
        # 1. Create the single-generation transition matrix for allele counts (0 to N)
        T_mat <- matrix(0, nrow = N + 1, ncol = N + 1)
        for (i in 0:N) {
            p <- i / N
            T_mat[i + 1, ] <- dbinom(0:N, size = N, prob = p)
        }

        # 2. Compute the multi-generation transition matrix using matrix exponentiation
        T_ngen <- expm::`%^%`(T_mat, n_gen)

        # 3. Define the VAF bin boundaries, identical to the simulation's logic
        bin_boundaries <- c(0, seq(0, 1, 1/k), 1)
        vaf_bins <- get_vaf_bins(k)

        # 4. Aggregate the allele count probabilities into the VAF bins
        A <- matrix(0, nrow = k + 2, ncol = k + 2)
        
        for (i in 1:(k + 2)) {
            # Determine the representative starting allele count for the i-th VAF bin,
            # using the bin midpoint, which mirrors the simulation's setup.
            start_vaf <- vaf_bins[i]
            start_allele_count <- round(start_vaf * N)
            start_allele_count <- max(0, min(N, start_allele_count))
            
            # Get the probability distribution of allele counts after n_gen generations
            prob_dist_after_ngen <- T_ngen[start_allele_count + 1, ]
            
            for (j in 1:(k + 2)) {
                # This binning logic now exactly replicates the original simulation function.
                xstart <- floor(bin_boundaries[j] * N) + 1
                xend <- floor(bin_boundaries[j+1] * N)
                xstart <- min(xstart, xend)

                if (j == k + 1) {
                    xend <- xend - 1
                }

                if (xstart > xend) {
                    A[i, j] <- 0
                } else {
                    indices <- xstart:xend
                    # Sum probabilities for allele counts in the current bin
                    A[i, j] <- sum(prob_dist_after_ngen[indices + 1])
                }
            }
        }
    }

    # 5. Apply the mutation rate 'eps' to the boundary conditions
    A[1, ] <- c(1 - eps, rep(eps / (k + 1), k + 1))
    A[nrow(A), ] <- rev(A[1, ])
    
    # Set row and column names for clarity
    colnames(A) <- vaf_bins
    rownames(A) <- vaf_bins

    return(A)
}

#' Wrapper function to interpolate non-integer generations
#' @param k number of VAF bins
#' @param eps error rate
#' @param N population size
#' @param n_gen number of generations
#' @return transition matrix
#' @export
get_transition_mat_wf_hmm_wrapper = function(k, eps, N, n_gen) {
    
    # Check if n_gen is an integer
    if (n_gen == round(n_gen)) {
        # If integer, use the original function directly
        return(get_transition_mat_wf_hmm(k = k, eps = eps, N = N, n_gen = n_gen))
    }
    
    # If non-integer, interpolate between two nearest integer generations
    n_gen_lower = floor(n_gen)
    n_gen_upper = ceiling(n_gen)
    
    # Get transition matrices for the two nearest integer generations
    A_lower = get_transition_mat_wf_hmm(k = k, eps = eps, N = N, n_gen = n_gen_lower)
    A_upper = get_transition_mat_wf_hmm(k = k, eps = eps, N = N, n_gen = n_gen_upper)
    
    # Calculate interpolation weight
    weight = n_gen - n_gen_lower
    
    # Linear interpolation between the two matrices
    A_interpolated = (1 - weight) * A_lower + weight * A_upper
    
    return(A_interpolated)
}

#' Convert a long format mutation data frame to a matrix
#' @param mut_dat_long A data frame with columns: variant, cell, variable
#' @param variable The variable to convert to a matrix
#' @return A matrix of values
#' @export
long_to_mat = function(mut_dat_long, variable) {
    as.data.table(mut_dat_long) %>%
        data.table::dcast(variant ~ cell, value.var = variable, fill = 0) %>%
        tibble::column_to_rownames("variant") %>%
        as.matrix()
}

#' Make a rooted NJ tree
#' @param vmat A matrix of cell-by-variable values
#' @param dist_method The distance method to use
#' @return A phylo object
#' @export
make_rooted_nj = function(vmat, dist_method = 'manhattan') {
    vmat = cbind(vmat, outgroup = 0)
    dist_mat = vmat %>% as.matrix %>% t %>% dist(method = dist_method)
    nj_tree = ape::nj(dist_mat) %>% 
        ape::root(outgroup = 'outgroup') %>% drop.tip('outgroup') 
    return(nj_tree)
}

modify_A = function(A, eps) {
    A[1,] = c(1-eps, rep(eps/(ncol(A)-1), ncol(A)-1))
    A[nrow(A),] = rev(A[1,])
    return(A)
}

scale_eps <- function(eps, f, mu = 1) {
  logit <- function(p) log(p / (1 - p))
  logistic <- function(x) 1 / (1 + exp(-x))
  logistic(logit(eps) + mu * f)
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


copy_crf = function(crf) {
    crf_copy <- rlang::env_clone(crf)
    attributes(crf_copy) <- attributes(crf)
    return(crf_copy)
}

#' @export 
decode_tree = function(
    tn, A, liks, post_max = FALSE, store_bels = FALSE, store_crfs = FALSE, debug = FALSE,
    score_only = FALSE
) {
    
    k = ncol(A)
    vafs = as.numeric(colnames(A))
    # convert tree to CRF
    tn$node.label = NULL
    Gn = as.igraph(tn)
    
    gtree = as_tbl_graph(Gn)
    root_node = gtree %>% filter(node_is_root()) %>% pull(name) %>% as.character

    adj_n = as_adjacency_matrix(Gn)
    crf = make.crf(adj_n, k)

    # add edge potentials
    flip_dict = Gn %>% 
        as_edgelist(names = F) %>% 
        apply(1, function(x){
            setNames(is.unsorted(x), paste0(sort(x), collapse = ','))
        }, simplify = F) %>%
        unlist

    for (i in 1:crf$n.edges) {
        
        epair = paste0(sort(crf$edges[i,]), collapse = ',')
        flip = flip_dict[[epair]]
        
        if (flip) {
            crf$edge.pot[[i]] = t(A)
        } else {
            crf$edge.pot[[i]] = A
        }
        
    }

    # add note potentials
    vnames = names(V(Gn))
    crf$node.labels = vnames
    rownames(crf$node.pot) = vnames

    logZ = c()
    ebels = list()
    nbels = list()
    crfs = list()

    for (mut in names(liks)) {

        crf$node.pot[colnames(liks[[mut]]),] = t(liks[[mut]])
        crf$node.pot[!vnames %in% colnames(liks[[mut]]),] = 1/k
        crf$node.pot[root_node,] = c(1, rep(0, k-1))
        
        # decoding
        res_mar = infer.tree(crf)
        
        if (!score_only) {
            # append posterior mean to graph tree
            p_dat = res_mar$node.bel %*% diag(vafs) %>% rowSums
            p_dat = p_dat %>% data.frame(vnames, .) %>%
                setNames(c('name', paste0('p_', mut)))
            
            gtree = gtree %>% activate(nodes) %>%
                select(-any_of(c(paste0('p_', mut)))) %>% 
                left_join(p_dat, by = join_by(name))

            if (post_max) {
                res_max = decode.tree(crf)
                z_dat = data.frame(vnames, vafs[res_max]) %>% setNames(c('name', paste0('z_', mut))) 

                gtree = gtree %>% activate(nodes) %>%
                    select(-any_of(c(paste0('z_', mut)))) %>% 
                    left_join(z_dat, by = join_by(name))
            } else {
                res_max = NULL
            }
        }

        logZ = c(logZ, res_mar$logZ)

        if (store_bels) {
            ebels[[mut]] = res_mar$edge.bel
            nbels[[mut]] = res_mar$node.bel
        }

        if (store_crfs) {
            crf_copy <- rlang::env_clone(crf)
            attributes(crf_copy) <- attributes(crf)
            crfs[[mut]] = crf_copy
        }
    }

    logZ = setNames(logZ, names(liks))
    gtree$logZ = logZ

    if (debug) {
        return(list('gtree' = gtree, 'crfs' = crfs, 'Gn' = Gn, 
        'res_mar' = res_mar, 'res_max' = res_max, 'ebels' = ebels, 'nbels' = nbels))
    }

    return(gtree)
}

# allow A to vary across mutations
#' @export 
decode_tree_brl = function(
    tn, As, liks, root_eps = 1e-2, post_max = FALSE, store_bels = FALSE, store_crfs = FALSE, debug = FALSE,
    score_only = FALSE
) {

    if (!inherits(tn, "igraph")) {
        Gn = as.igraph(tn)
        E(Gn)$length = tn$edge.length
    } else {
        Gn = tn
        if (is.null(E(Gn)$length)) {
            stop('No edge length in tree')
        }
    }
    
    k = ncol(As[[1]])
    vafs = as.numeric(colnames(As[[1]]))
    
    gtree = as_tbl_graph(Gn)
    root_node = gtree %>% filter(node_is_root()) %>% pull(name) %>% as.character

    adj_n = as_adjacency_matrix(Gn)
    crf = make.crf(adj_n, k)

    # add edge potentials
    elist = Gn %>% as_edgelist(names = F)

    flip_df = data.frame(
            from = pmin(elist[,1], elist[,2]),
            to = pmax(elist[,1], elist[,2]),
            flip = elist[,1] > elist[,2],
            length = E(Gn)$length
        ) %>%
        arrange(from, to)

    crf$edge.pot = lapply(
        1:nrow(flip_df),
        function(i){

            flip = flip_df[i,]$flip
            n_gen = flip_df[i,]$length

            A = interpolate_matrices(As, n_gen, min_index = 1e-5)

            if (flip) {
                t(A)
            } else {
                A
            }

    })

    # add note potentials
    vnames = names(V(Gn))
    crf$node.labels = vnames
    rownames(crf$node.pot) = vnames

    logZ = c()
    ebels = list()
    nbels = list()
    crfs = list()

    for (mut in names(liks)) {

        crf$node.pot[colnames(liks[[mut]]),] = t(liks[[mut]])
        crf$node.pot[!vnames %in% colnames(liks[[mut]]),] = 1/k
        # root node doesn't have to clean
        crf$node.pot[root_node,] = c(1-root_eps, rep(root_eps/(k-1), k-1))
        
        # decoding
        res_mar = infer.tree(crf)
        
        if (!score_only) {
            # append posterior mean to graph tree
            p_dat = res_mar$node.bel %*% diag(vafs) %>% rowSums
            p_dat = p_dat %>% data.frame(vnames, .) %>%
                setNames(c('name', paste0('p_', mut)))
            
            gtree = gtree %>% activate(nodes) %>%
                select(-any_of(c(paste0('p_', mut)))) %>% 
                left_join(p_dat, by = join_by(name))

            if (post_max) {
                res_max = decode.tree(crf)
                z_dat = data.frame(vnames, vafs[res_max]) %>% setNames(c('name', paste0('z_', mut))) 

                gtree = gtree %>% activate(nodes) %>%
                    select(-any_of(c(paste0('z_', mut)))) %>% 
                    left_join(z_dat, by = join_by(name))
            } else {
                res_max = NULL
            }
        }

        logZ = c(logZ, res_mar$logZ)

        if (store_bels) {
            ebels[[mut]] = res_mar$edge.bel
            nbels[[mut]] = res_mar$node.bel
        }

        if (store_crfs) {
            crf_copy <- rlang::env_clone(crf)
            attributes(crf_copy) <- attributes(crf)
            crfs[[mut]] = crf_copy
        }
    }

    logZ = setNames(logZ, names(liks))
    gtree$logZ = logZ

    if (debug) {
        return(list('gtree' = gtree, 'crfs' = crfs, 'Gn' = Gn, 
        'res_mar' = res_mar))
    }

    return(gtree)
}

interpolate_matrices <- function(As, index, min_index = 1e-5) {

    index = max(min_index, index)

    # Get the list of provided indices
    provided_indices <- as.numeric(names(As))
    min_index <- min(provided_indices)
    max_index <- max(provided_indices)

    # If index is below the minimum, return the matrix at the minimum index
    if (index <= min_index) {
        return(As[[as.character(min_index)]])
    }

    # If index is above the maximum, return the matrix at the maximum index
    if (index >= max_index) {
        return(As[[as.character(max_index)]])
    }


    # If the index exactly matches one of the provided indices, return the corresponding matrix directly
    if (index %in% provided_indices) return(As[[as.character(index)]])

    # Find the closest lower and upper indices for interpolation
    lower_index <- max(provided_indices[provided_indices < index])
    upper_index <- min(provided_indices[provided_indices > index])

    # Calculate the weight for interpolation
    weight <- (index - lower_index) / (upper_index - lower_index)

    # Retrieve matrices for interpolation
    A_lower <- As[[as.character(lower_index)]]
    A_upper <- As[[as.character(upper_index)]]

    # Ensure the matrices are the same dimensions
    if (!all(dim(A_lower) == dim(A_upper))) stop("Matrices must have the same dimensions.")

    # Perform linear interpolation
    interpolated_matrix <- (1 - weight) * A_lower + weight * A_upper

    # Normalize each row to sum to 1
    interpolated_matrix <- interpolated_matrix / rowSums(interpolated_matrix)

    return(interpolated_matrix)
}

estimate_brl = function(phy, liks, As, height = NULL, init = 1) {

    if (length(init) == 1) {
        init = rep(init, length(phy$edge.length))
    }

    G = igraph::as.igraph(phy)
    
    fit = optim(
        fn = function(x) {

            phy$edge.length = x
            
            gtree = decode_tree_brl(
                tn = phy,
                As,
                liks, score_only = T)
            
            E(G)$weight = x
            dists = distances_to_root_tips(G)

            if (is.null(height)) {
                height = mean(dists)
            }

            mse = mean((dists - height)^2)
            
            -sum(gtree$logZ) + mse
            
        },
        method = 'L-BFGS-B',
        par = init,
        lower = 1,
        upper = 100
    )

    phy$edge.length = fit$par

    return(phy)
}

# Use internal branch representation for ultrametric tree
# see https://github.com/NickWilliamsSanger/rtreefit
estimate_brl_recode = function(phy, As, liks, htotal, init = 0.1) {

    if (length(init) == 1) {
        init = rep(init, phy$Nnode)
    }

    phy = reorder(phy)

    G = igraph::as.igraph(phy)
    
    fit = optim(
        fn = function(x) {

            G = x_to_length(G, x, htotal = htotal)
            
            gtree = decode_tree_brl(
                tn = G,
                As,
                liks, score_only = T)
        
            -sum(gtree$logZ)
            
        },
        method = 'L-BFGS-B',
        par = init,
        lower = 0,
        upper = 1
    )

    G = x_to_length(G, fit$par, htotal = htotal)

    phy$edge.length = E(G)$length

    return(phy)
}


# Use internal branch representation for ultrametric tree
# see https://github.com/NickWilliamsSanger/rtreefit
#' @export 
estimate_brl_par = function(tree, As, liks, htotal, x_init = 0.1, change_h = FALSE) {

    if (!'tbl_graph' %in% class(tree)) {
        tree = reorder(tree)
        G = igraph::as.igraph(tree)
    } else {
        G = tree
    }

    n_nodes = sum(degree(G, mode = "all") > 1)
    
    if (length(x_init) == 1) {
        x_init = rep(x_init, n_nodes)
    }

    if (change_h) {

        fn = function(x) {
            
            G = x_to_length(G, x[-1], htotal = x[1])
            
            gtree = decode_tree_brl(
                tn = G,
                As,
                liks, score_only = TRUE)
        
            -sum(gtree$logZ)
            
        }

        fit = optimParallel(
            fn = fn,
            method = 'L-BFGS-B',
            par = c(htotal, x_init),
            lower = c(1, rep(0, n_nodes)),
            upper = c(2000, rep(1, n_nodes))
        )

        G = x_to_length(G, fit$par[-1], htotal = fit$par[1])
        
    } else {

        fn = function(x) {

            G = x_to_length(G, x, htotal = htotal)
            
            gtree = decode_tree_brl(
                tn = G,
                As,
                liks, score_only = TRUE)
        
            -sum(gtree$logZ)
            
        }

        fit = optimParallel(
                fn = fn,
                method = 'L-BFGS-B',
                par = x_init,
                lower = 0,
                upper = 1
            )

        G = x_to_length(G, fit$par, htotal = htotal)

    }

    return(G)
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
                phylist = lapply(
                    elist,
                    function(edges){
                    phy_new = attach_edges(phy, edges)
                    return(phy_new)
                })

                class(phylist) = 'multiPhylo'

                return(phylist)
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
    move_type = 'NNI_multi'
) {

    # Rearrangements have to be unrooted
    if (move_type == 'NNI') {
        propose_tree = TreeSearch::NNI
    } else if (move_type == 'NNI_multi') {
        propose_tree = rNNI_multi
    } else if (move_type == 'NNT') {
        propose_tree = rNNT
    } else {
        stop('Invalid move type')
    }

    res = mclapply(
        1:nchains,
        mc.cores = ncores,
        function(s) {

            set.seed(s)

            # need to change this to pre-allocation
            phy_list = vector("list", max_iter + 1)
            l_0 = sum(gtree_init$logZ)
            phy_list[[1]] = as.phylo(gtree_init)

            for (i in 1:max_iter) {
                
                phy_new = propose_tree(phy_list[[i]])
                gtree_new = decode_tree(phy_new, A, liks, score_only = TRUE)
                
                l_1 = sum(gtree_new$logZ)
                probab = exp(l_1 - l_0)

                if (runif(1) < probab) {
                    phy_list[[i+1]] = phy_new
                    l_0 = l_1
                } else {
                    phy_list[[i+1]] = phy_list[[i]]
                }
            }

            if (trace &!is.null(outfile)) {
                outdir = dirname(outfile)
                fname = basename(outfile)
                saveRDS(phy_list, glue('{outdir}/chain{s}_{fname}'))
            }

            return(phy_list)
        }
    )

    return(res)
}

# to fix: apparently the init tree has to be rooted otherwise to_phylo_reoder won't work. 
#' @export
run_tree_mcmc = function(
    phy_init, logP_list, logA_vec, max_iter = 100, nchains = 1, ncores = 1, outfile = NULL, resume = FALSE
) {

    chains = 1:nchains

    if (!is.null(outfile)) {
        outdir = dirname(outfile)
        fname = basename(outfile)

        if (!dir.exists(outdir)) {
            dir.create(outdir, recursive = TRUE)
        }
        
        if (resume) {
            done = sapply(chains, function(s) {file.exists(glue('{outdir}/chain{s}_{fname}'))})
            chains = chains[!done]
        }
    }

    message('Running MCMC with ', length(chains), ' chains')



    res = mclapply(
        chains,
        mc.cores = ncores,
        function(s) {

            elist = tree_mcmc_cpp(phy_init$edge, logP_list, logA_vec, max_iter = max_iter, seed = s)

            tree_list = lapply(
                elist,
                function(edges){
                    attach_edges(phy_init, edges)
            })

            class(tree_list) = 'multiPhylo'

            if (!is.null(outfile)) {
                saveRDS(tree_list, glue('{outdir}/chain{s}_{fname}'))
            }

            return(tree_list)
        }
    )

    if (resume & !is.null(outfile)) {

        chains_prev = setdiff(1:nchains, chains)

        res_prev = mclapply(
            chains_prev,
            mc.cores = ncores,
            function(s) {
                readRDS(glue('{outdir}/chain{s}_{fname}'))
            }
        )
        res = c(res_prev, res)

    }

    saveRDS(res, outfile)

    return(res)
}

#' @export
collect_chains = function(res, burnin = 0, max_iter = Inf) {

    if (max_iter < burnin) {
        stop('Max iter needs to be greater than burnin')
    }

    res = res[res %>% sapply(length) > 0]

    mcmc_trees = res %>% lapply(function(trees){
            trees = trees[(burnin+1):min(length(trees), max_iter)]
            lapply(trees, function(tree){
                tree$edge.length = NULL
                tree$node.label = NULL
                tree$nodes = NULL
                tree$edge = TreeTools::RenumberTree(tree$edge[,1], tree$edge[,2])
                return(tree)
            })
        }) %>% unlist(recursive = F)

    class(mcmc_trees) = 'multiPhylo'

    return(mcmc_trees)
    
}

#' @export
add_conf = function(gtree, phylist) {

    if (!'multiPhylo' %in% class(phylist) & is.list(phylist)) {
        class(phylist) = 'multiPhylo'
    }

    if (!'tbl_graph' %in% class(gtree)) {
        gtree$node.label = NULL
        gtree = as_tbl_graph(gtree)
    }
    
    conf_dict = to_phylo_reorder(gtree) %>%
        add_clade_freq(phylist) %>%
        parse_conf()

    gtree = gtree %>% activate(nodes) %>% mutate(conf = conf_dict[name])

    return(gtree)
    
}

# reorder and preserves node labels
#' @export
to_phylo_reorder = function(graph) {
    df <- igraph::as_data_frame(graph)
    node_counts <- table(c(df$to, df$from))
    tips <- names(node_counts)[node_counts == 1]
    nodes <- names(node_counts)[node_counts > 1]
    attr <- igraph::vertex_attr(graph)
    tipn <- 1:length(tips)
    names(tipn) <- tips
    noden <- (length(tips) + 1):(length(tips) + length(nodes))
    names(noden) <- nodes
    renumber <- c(tipn, noden)
    df$from <- as.numeric(renumber[df$from])
    df$to <- as.numeric(renumber[df$to])
    phylo <- list()
    phylo$edge <- matrix(cbind(df$from, df$to), ncol = 2)
    phylo$edge.length <- as.numeric(df$length)
    phylo$tip.label <- tips
    phylo$node.label <- nodes
    phylo$Nnode <- length(nodes)
    class(phylo) <- "phylo"
    nnodes <- length(renumber)
    phylo$nodes <- lapply(1:nnodes, function(x) {
        n <- list()
        n$id <- names(renumber[renumber == x])
        n
    })
    phylo = reorder(phylo)
    if (length(phylo$edge.length) == 0) {
        phylo$edge.length = NULL
    }
    return(phylo)
}

#' @export
add_clade_freq = function(phy, phys, rooted = TRUE) {
    freqs = prop.clades(phy, phys, rooted = rooted)/length(phys)
    phy$node.label = freqs
    return(phy)
}

#' @export
parse_conf = function(phy) {
    nodes = unname(unlist(phy$nodes))
    nodes = nodes[!nodes %in% phy$tip.label]
    conf_dict = setNames(phy$node.label, nodes)
    return(conf_dict)
}

#' @export
phylo_to_gtree = function(phy) {
        
    tip_nodes <- data.frame(
        name = phy$tip.label
    )
    
    internal_nodes <- data.frame(
        name = paste0("Node", seq_len(phy$Nnode))
    )

    if (!is.null(phy$node.label)) {
        internal_nodes$label = phy$node.label
    }
    
    nodes <- bind_rows(tip_nodes, internal_nodes)
    edges <- data.frame(phy$edge)
    colnames(edges) <- c("from", "to")
    
    gtree <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
    
    return(gtree)
}

#' @export
trim_tree = function(tree, conf) {
    tree$node.label[is.na(tree$node.label)] = 0
    tree = TreeTools::CollapseNode(tree, which(tree$node.label < conf) + length(tree$tip.label))
    tree = TreeTools::Renumber(tree)
    return(tree)
}

#' @export
get_consensus = function(phylist, p = 0.5, check_labels = FALSE, rooted = TRUE, conf = FALSE) {
    phy_cons = ape::consensus(phylist, p = p, rooted = rooted, check.labels = check_labels)
    gtree_cons = phylo_to_gtree(phy_cons) %>% rename(conf = label)
    return(gtree_cons)
}

#' @export
getConfidentClades <- function(pp, p = 0.9, max_size = Inf, labels = TRUE, singletons = TRUE) {
    # extract clade sizes and confidences
    tip_labels   <- attr(pp, "labels")
    sizes <- vapply(pp, length, integer(1))
    nums  <- attr(pp, "number")
    nt    <- nums[1]
    if (is.null(nums)) {
        stop("prop.part object must have 'number' attribute")
    }
    confs <- nums / nt

    # exclude the trivial partition (all cells) and too‐large clades
    trivial_idx   <- which(sizes == max(sizes))
    candidate_idx <- setdiff(seq_along(pp), trivial_idx)
    candidate_idx <- candidate_idx[sizes[candidate_idx] <= max_size]

    # 1) iterative selection of confident, disjoint clades
    selected_idx <- list()
    remaining   <- candidate_idx
    repeat {
        # rank remaining by descending size
        ord  <- remaining[order(sizes[remaining], decreasing = TRUE)]
        # pick first with conf > p
        keep <- ord[confs[ord] > p]
        if (length(keep) == 0) break
        pick <- keep[1]
        selected_idx[[length(selected_idx) + 1]] <- pp[[pick]]
        # remove any clades overlapping this pick
        pick_cells <- pp[[pick]]
        remaining  <- remaining[
            vapply(remaining, function(i) {
                length(intersect(pp[[i]], pick_cells)) == 0
            }, logical(1))
        ]
    }

    # 2) identify all tips and add singletons for any unassigned tips
    if (singletons) {
        all_tips_idx <- seq_along(tip_labels)
        assigned     <- unique(unlist(selected_idx))
        unassigned   <- setdiff(all_tips_idx, assigned)
        # add each unassigned tip as a singleton clade
        selected_idx <- c(
            selected_idx,
            lapply(unassigned, function(i) i)
        )
    }

    # 3) optionally map indices to labels
    if (labels) {
        selected <- lapply(selected_idx, function(cl) tip_labels[cl])
    } else {
        selected <- selected_idx
    }

    return(selected)
}


# write a function to estimate parameters using MCMC with fmcmc package
#' @param tree_fit phylogenetic tree
#' @param amat adjacency matrix
#' @param dmat distance matrix
#' @param initial_params initial parameter values (n_gen, log_eps, log_err)
#' @param lower_bounds lower bounds for parameters in transformed space
#' @param upper_bounds upper bounds for parameters in transformed space
#' @param nsteps number of MCMC steps
#' @param nchains number of MCMC chains
#' @param n_pop population size for likelihood computation
#' @param k number of clusters for likelihood computation
#' @param ncores number of cores to use
#' @param burnin number of MCMC steps to discard as burnin
#' @param outfile file to save MCMC result
#' @return MCMC result object from fmcmc
#' @export
fit_params_mcmc = function(
    tree_fit, amat, dmat,
    initial_params = c('n_gen' = 100, 'log_eps' = log(1e-3), 'log_err' = log(1e-3)),
    lower_bounds = c('n_gen' = 1, 'log_eps' = log(1e-18), 'log_err' = log(1e-18)),
    upper_bounds = c('n_gen' = 1000, 'log_eps' = log(0.2), 'log_err' = log(0.2)),
    nsteps = 500, nchains = 1, outfile = NULL,
    n_pop = 600, k = 20, ncores = 1, burnin = 100) {
    
    # Ensure tree is properly formatted
    tree_fit = reorder_phylo(tree_fit)
    tree_fit$edge.length = NULL
    
    # All three parameters will be estimated
    param_names = c('n_gen', 'log_eps', 'log_err')
    
    message(glue("Estimating all 3 parameters using MCMC: {paste(param_names, collapse=', ')}"))
    
    # Create log-likelihood function for MCMC
    # This function takes parameters in transformed space and converts them back
    log_likelihood = function(p, tree_fit., amat., dmat., n_pop., k.) {
        # Convert parameters back to original scale
        params_orig = p
        
        # Convert log-scale parameters back to original scale
        params_orig['log_eps'] = exp(p['log_eps'])
        names(params_orig)[names(params_orig) == 'log_eps'] = 'eps'
        params_orig['log_err'] = exp(p['log_err'])
        names(params_orig)[names(params_orig) == 'log_err'] = 'err'
        
        # Compute likelihood
        ll = get_param_lik_cpp(tree_fit., amat., dmat., 
            n_gen = params_orig['n_gen'],
            eps = params_orig['eps'], 
            err = params_orig['err'],
            n_pop = n_pop., k = k.)
        
        return(ll)
    }
    
    kernel = fmcmc::kernel_normal_reflective(
            scale = c(5, 0.1, 0.1),
            scheme = "joint",
            lb = lower_bounds,
            ub = upper_bounds
        )
    
    # Run MCMC
    message("Starting MCMC sampling...")

    res_df_all = mclapply(
        1:nchains,
        mc.cores = ncores,
        function(i) {
            set.seed(i)
                    
            mcmc_result = fmcmc::MCMC(
                log_likelihood,
                initial = initial_params,
                nsteps = nsteps,
                nchains = 1,
                kernel = kernel,
                tree_fit. = tree_fit,
                amat. = amat,
                dmat. = dmat,
                n_pop. = n_pop,
                k. = k
            )

            res_df = data.frame(mcmc_result) %>% 
                tibble::rowid_to_column('iter') %>%
                mutate(
                    err = exp(log_err), 
                    eps = exp(log_eps),
                    chain = i
                ) %>%
                reshape2::melt(id.var = c('iter', 'chain'))

            return(res_df)
        }
    ) %>% bind_rows()

    if (!is.null(outfile)) {
        fwrite(res_df_all, outfile)
    }

    params_est = res_df_all %>%
        filter(iter > burnin) %>%
        group_by(variable) %>%
        summarise(est = mean(value), .groups = 'drop') %>%
        {setNames(.$est, .$variable)}

    return(params_est)
}