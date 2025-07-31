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

#' @export
get_transition_mat_wf_hmm <- function(k, eps = 0.01, N = 100, n_gen = 100) {
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

nnin_multi <- function(tree, n, d = Inf) {
  # Remove any existing "order" attribute so reorder() works cleanly
  attr(tree, "order") <- NULL
  
  # Prepare variables
  edge         <- matrix(tree$edge, ncol = 2)
  parent       <- edge[, 1]
  child        <- edge[, 2]
  leaf_cutoff  <- min(parent) - 1  # tips are 1:leaf_cutoff
  
  # Find the n-th internal edge (child > leaf_cutoff means "child is an internal node")
  ind <- which(child > leaf_cutoff)[n]
  if (is.na(ind)) return(NULL)
  
  # p1 → p2 is the selected internal edge
  p1 <- parent[ind]
  p2 <- child[ind]
  
  # All children of p1 (excluding p2)
  c1_inds    <- which(parent == p1)
  siblings1  <- child[c1_inds]
  x_list     <- setdiff(siblings1, p2)
  
  # All children of p2 (excluding p1, if p1 appears in that list)
  c2_inds <- which(parent == p2)
  y_list  <- child[c2_inds]
  
  # Randomly sample up to d children from each side
  if (length(x_list) > d) {
    x_list <- sample(x_list, d)
  }
  if (length(y_list) > d) {
    y_list <- sample(y_list, d)
  }
  
  # Container for resulting trees
  result <- list()
  
  # For each sampled pair (x, y), swap x under p2 and y under p1
  for (x in x_list) {
    for (y in y_list) {
      tr_new <- tree
      
      # Locate the edge-rows where child == x and child == y
      idx_x <- which(tr_new$edge[, 2] == x)
      idx_y <- which(tr_new$edge[, 2] == y)
      
      # Swap parents: detach x from p1 and attach to p2; detach y from p2 and attach to p1
      tr_new$edge[idx_x, 1] <- p2
      tr_new$edge[idx_y, 1] <- p1
      
      # If there are edge lengths, swap those too so lengths stay consistent
      if (!is.null(tree$edge.length)) {
        tmp <- tr_new$edge.length[idx_x]
        tr_new$edge.length[idx_x] <- tr_new$edge.length[idx_y]
        tr_new$edge.length[idx_y] <- tmp
      }
      
      # Re-order the tree in postorder so ape will maintain a valid phylo object
      tr_new <- reorder(tr_new, "postorder")
      
      # Add to result list
      result[[length(result) + 1]] <- tr_new
    }
  }
  
  # Return as multiPhylo
  class(result) <- "multiPhylo"
  return(result)
}

nni_multi <- function(tree, d = Inf) {
  # 1) Identify internal edges by finding rows where 'child > k'
  edge   <- matrix(tree$edge, ncol = 2)
  parent <- edge[, 1]
  child  <- edge[, 2]
  k      <- min(parent) - 1
  
  # internal_edges is the vector of row‐indices for edges whose 'child' is an internal node
  internal_edges <- which(child > k)
  if (length(internal_edges) == 0) {
    return(list())  # no internal edges → no NNI moves possible
  }
  
  # 2) For each internal edge (indexed by its position in the 'child>k' list),
  #    call nnin_multi(tree, n) where n = rank among internal edges
  results <- list()
  for (i in seq_along(internal_edges)) {
    # i is the rank of this internal edge among all internal edges
    candidate_trees <- nnin_multi(tree, i, d)
    if (!is.null(candidate_trees) && length(candidate_trees) > 0) {
      # Append all returned trees to our master list
      for (tr in candidate_trees) {
        results[[ length(results) + 1 ]] <- tr
      }
    }
  }

  class(results) = 'multiPhylo'
  
  return(results)
}

rNNI_multi <- function(tree, d = Inf) {
  # Identify “k” so that any child > k is an internal node
  edge   <- matrix(tree$edge, ncol = 2)
  parent <- edge[,1]
  child  <- edge[,2]
  k      <- min(parent) - 1
  
  # Find all internal edges (rows where child > k)
  internal_edges <- which(child > k)
  if (length(internal_edges) == 0) {
    return(NULL)  # no internal edges → no NNI moves
  }
  
  # Randomly pick one of those internal edges by its rank n
  # (nnin_multi expects “n” to be the index among child>k)
  n <- sample(seq_along(internal_edges), 1)
  
  # Get the list of possible swaps for that edge
  moves <- nnin_multi(tree, n, d)
  if (is.null(moves) || length(moves) == 0) {
    return(NULL)  # edge had no valid multifurcation‐NNI moves
  }
  
  # Randomly choose one of the returned trees
  selected <- moves[[ sample(seq_along(moves), 1) ]]
  return(selected)
}


nntn <- function(tree, n) {
  # ——— remove any existing “order” attribute, so reorder() won’t choke ———
  attr(tree, "order") <- NULL

  # Extract parent/child arrays:
  edges       <- matrix(tree$edge, ncol = 2)
  parent      <- edges[,1]
  child       <- edges[,2]
  leaf_cutoff <- min(parent) - 1    # nodes 1:leaf_cutoff are tips

  # All “internal edges” are those where child > leaf_cutoff
  internal_edges <- which(child > leaf_cutoff)
  if (n < 1 || n > length(internal_edges)) {
    return(list())  # invalid index → no neighbors
  }
  ind <- internal_edges[n]

  # Let p1 -> p2 be the chosen internal edge
  p1 <- parent[ind]
  p2 <- child[ind]

  # Gather p1’s children (including p2) and p2’s children
  c1_inds   <- which(parent == p1)
  children1 <- child[c1_inds]      # all nodes directly under p1
  x_list    <- setdiff(children1, p2)  # all siblings of p2 under p1

#   print(x_list)

  c2_inds   <- which(parent == p2)
  children2 <- child[c2_inds]      # all nodes directly under p2
  y_list    <- children2           # all children under p2

#   print(y_list)

  result <- list()

  # 1) Move “x” from p1 → p2 (downwards), but only if p1 has ≥ 3 children
  #    (so that after removing x, p1 still has ≥2 children)
  if (length(children1) >= 3) {
    for (x in x_list) {
      tr_new <- tree

      # re‐attach x under p2
      idx_x <- which(tr_new$edge[,2] == x)
      tr_new$edge[idx_x, 1] <- p2

      # reorder to get a valid “phylo” object again
      tr_new <- reorder(tr_new, "postorder")
      result[[ length(result) + 1 ]] <- tr_new
    }
  }

  # 2) Move “y” from p2 → p1 (upwards), but only if p2 has ≥ 3 children
  #    (so that after removing y, p2 still has ≥2 children)
  if (length(children2) >= 3) {
    for (y in y_list) {
      tr_new <- tree

      # re‐attach y under p1
      idx_y <- which(tr_new$edge[,2] == y)
      tr_new$edge[idx_y, 1] <- p1

      tr_new <- reorder(tr_new, "postorder")
      result[[ length(result) + 1 ]] <- tr_new
    }
  }

  class(result) <- "multiPhylo"
  return(result)
}


nnt <- function(tree) {
  # Remove any existing “order” attribute
  attr(tree, "order") <- NULL

  # Identify all internal edges
  edges       <- matrix(tree$edge, ncol = 2)
  parent      <- edges[,1]
  child       <- edges[,2]
  leaf_cutoff <- min(parent) - 1
  internal_edges <- which(child > leaf_cutoff)

  if (length(internal_edges) == 0) {
    return(list())  # no internal edges → no moves
  }

  all_neighbors <- list()
  for (i in seq_along(internal_edges)) {
    nb_i <- nntn(tree, i)
    if (length(nb_i) > 0) {
      for (tr in nb_i) {
        all_neighbors[[ length(all_neighbors) + 1 ]] <- tr
      }
    }
  }

  class(all_neighbors) <- "multiPhylo"
  return(all_neighbors)
}

rNNT <- function(tree) {
  all_nb <- nnt(tree)  
  choice <- sample(seq_along(all_nb), 1)
  return(all_nb[[choice]])
}

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

find_far_nodes = function(phy) {
    # assume phy is a rooted ape::phylo
    root <- Ntip(phy) + 1
    
    # 1) immediate internal children of the root
    imm_int <- phy$edge[phy$edge[,1] == root, 2]
    imm_int <- imm_int[imm_int > Ntip(phy)]
    
    # 2) all internal nodes excluding the root
    all_int <- (Ntip(phy)+1):(Ntip(phy) + phy$Nnode)
    all_int = all_int[all_int != root]
    
    # 3) those at distance > 1 from root
    far_int <- setdiff(all_int, imm_int)

    return(far_int)

}

merge_trees = function(subtrees, collapse = FALSE) {
    phy_merge = subtrees %>%
        lapply(function(tr){TreeTools::AddTip(tr, label = 'out', where = 0)}) %>%
        Reduce(bind.tree, .) %>%
        drop.tip('out')
    if (collapse) {
        phy_merge = phy_merge %>% TreeTools::CollapseNode(find_far_nodes(.))
    }
    return(phy_merge)
}

one_node_depth2 <- function(tree, exclude_root = FALSE) {
    bifur  <- castor::is_bifurcating(tree)
    Ntips  <- length(tree$tip.label)
    root   <- Ntips + 1L

    # find all edges that end in a tip
    tipnodes <- which(tree$edge[,2] %in% seq_len(Ntips))
    edges    <- tree$edge[tipnodes, , drop=FALSE]

    df <- data.frame(
        node   = edges[,1],
        tip    = edges[,2],
        pat    = tree$edge.length[tipnodes]
    )

    if (exclude_root) {
        df <- df %>% dplyr::filter(node != root)
    }

    # now group by internal node and keep only those with ≥2 tips
    df2 <- df %>%
        dplyr::group_by(node) %>%
        dplyr::mutate(clades = list(tip)) %>%
        dplyr::rowwise() %>%
        dplyr::filter(length(clades) >= 2)

    # build the tip‐tip adjacency
    indices <- df2 %>%
        dplyr::select(-tip, -pat) %>%
        dplyr::distinct() %>%
        dplyr::pull(clades) %>%
    {
        if (!bifur) {
            lapply(., function(x) t(utils::combn(x, 2))) %>%
                do.call(rbind, .)
        } else {
            do.call(rbind, .)
        }
    }

    mat <- Matrix::sparseMatrix(
        i = c(indices[,1], indices[,2]),
        j = c(indices[,2], indices[,1]),
        dims = c(Ntips, Ntips)
    )

    mp <- mean(df2$pat) * 2

    list(W = mat, mean.pat = mp, pats = df2$pat)
}

one_node_tree_dist2 = function(tree, norm = TRUE, exclude_root = FALSE) {
    w <- one_node_depth2(tree, exclude_root = exclude_root)$W
    if (norm == TRUE) {
        w <- PATH:::rowNorm(w)
    }
    w <- w/sum(w)
    return(w)
}

# apparently rooted and unrooted produce different number of internal nodes
extended_consensus <- function(trees = NULL, pp = NULL, p = 1, check.labels = TRUE, rooted = FALSE) {

  # Recursive helper to rebuild topology from clusters
  foo <- function(ic, node) {
    pool <- pp[[ic]]
    if (ic < m) {
      for (j in (ic + 1):m) {
        wh <- match(pp[[j]], pool)
        if (!any(is.na(wh))) {
          edge[pos, 1] <<- node
          pool       <- pool[-wh]
          edge[pos, 2] <<- nextnode <<- nextnode + 1L
          pos        <<- pos + 1L
          foo(j, nextnode)
        }
      }
    }
    # Attach any remaining tips directly
    size <- length(pool)
    if (size) {
      ind <- pos:(pos + size - 1)
      edge[ind, 1] <<- node
      edge[ind, 2] <<- pool
      pos <<- pos + size
    }
  }

  if (!is.null(trees)) {
    # Align tip labels
    if (!is.null(attr(trees, "TipLabel"))) {
        labels <- attr(trees, "TipLabel")
    } else {
        labels <- trees[[1]]$tip.label
        if (check.labels) trees <- .compressTipLabel(trees)
    }

    # Optionally unroot for unrooted consensus
    if (!rooted) trees <- root(trees, 1)

    ntree <- length(trees)

    # Extract all splits and their frequencies
    if (is.null(pp)) {
        pp <- prop.part(trees, check.labels = FALSE)
    }
  }

  labels <- attr(pp, "labels")
  ntree  <- attr(pp, "number")[1]

  if (!rooted) {
    pp <- postprocess.prop.part(pp, "SHORTwise")
    pp[[1]] <- seq_along(labels)
  }

  # Threshold by frequency
  bs  <- attr(pp, "number")
  sel <- bs >= p * ntree
  pp  <- pp[sel]; bs <- bs[sel]

  # Remove trivial single-tip clusters
  lens <- lengths(pp)
  drop <- which(lens == 1)
  if (length(drop)) {
    pp  <- pp[-drop]
    bs  <- bs[-drop]
    lens <- lens[-drop]
  }

  # Sort by cluster size
  ind <- order(lens, decreasing = TRUE)
  pp  <- pp[ind]; bs <- bs[ind]

  # Greedy compatibility filter
  keep_idx <- integer(0)
  for (i in seq_along(pp)) {
    cl <- pp[[i]]
    ok <- TRUE
    for (j in keep_idx) {
      sj <- pp[[j]]
      ints <- intersect(cl, sj)
      if (!(length(ints) == 0 || length(ints) == length(cl) || length(ints) == length(sj))) {
        ok <- FALSE; break
      }
    }
    if (ok) keep_idx <- c(keep_idx, i)
  }
  pp <- pp[keep_idx]; bs <- bs[keep_idx]
  m  <- length(pp)

  # Build edge matrix
  n    <- length(labels)
  edge <- matrix(0L, n + m - 1, 2)
  if (m == 1) {
    # star tree: all tips attach to single node
    edge[, 1] <- n + 1L
    edge[, 2] <- 1:n
  } else {
    nextnode <- n + 1L
    pos      <- 1L
    foo(1, nextnode)
  }

  # Assemble phylo object
  res <- structure(
    list(edge = edge, tip.label = labels, Nnode = m),
    class = "phylo"
  )
  res <- reorder(res)

  # Attach support values to nodes
  res$node.label <- bs / ntree
  res
}

