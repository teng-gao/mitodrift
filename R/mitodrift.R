#' @import dplyr
#' @import tidygraph
#' @import stringr
#' @import fmcmc
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
    tree_init = NULL, logP, logA, max_iter = 100, 
    outfile = NULL, resume = FALSE, ncores = 1, trace_interval = 5
) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    if (is.list(logA)) {
        score_tree_func = score_tree_bp_wrapper_multi
        nni_func = nni_cpp_parallel_multi
    } else {
        score_tree_func = score_tree_bp_wrapper2
        nni_func = nni_cpp_parallel_cached
    }

    if (resume & !is.null(outfile) & file.exists(outfile)) {
        message('Resuming from existing tree list')
        tree_list = safe_read_chain(outfile)
        tree_current = tree_list %>% .[[length(.)]]
        max_current = tree_current$logZ
        start_iter = length(tree_list)
    } else {
        if (is.null(tree_init) || is.null(logP) || is.null(logA)) {
            stop("tree_init, logP, and logA must be provided when resume is FALSE")
        }
        start_iter = 1
        tree_init = reorder_phylo(tree_init)
        tree_init$edge.length = NULL
        tree_current = tree_init
        max_current = sum(score_tree_func(tree_current$edge, logP, logA))
        tree_current$logZ = max_current
        tree_list = phytools::as.multiPhylo(tree_current)
    }

    runtime = c(0,0,0)

    if (start_iter > max_iter) {
        message(paste0("Already completed ", start_iter - 1, " iterations. No new iterations to run for max_iter = ", max_iter, "."))
        if (!is.null(outfile)) {
            qs2::qd_save(tree_list, outfile)
        }
        return(tree_list)
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

        tree_list = tree_list %>% c(phytools::as.multiPhylo(tree_current))

        if (!is.null(outfile)) {
            if (i == 1 | i %% trace_interval == 0) {
                qs2::qd_save(tree_list, outfile)
            }
        }

        runtime = proc.time() - ptm
        
    }
    
    if (!is.null(outfile)) {
        qs2::qd_save(tree_list, outfile)
    }

    return(tree_list)
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
                    qs2::qd_save(tree_list, outfile)
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
                                dbinom(x = x$a, size = x$d, prob = pmin(v + eps, 1 - eps), log = log)
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
	# Precompute clamped probabilities for all VAF bins once
	p <- pmin(vafs + eps, 1 - eps)
	K <- length(p)
	variants <- rownames(amat)
	ncells <- ncol(amat)
	cell_names <- colnames(amat)

	# Use at most one core per variant
	nc <- min(ncores, length(variants))

	liks <- mclapply(
		variants,
		mc.cores = nc,
		function(v) {
			# Extract counts for this variant (length = #cells)
			x <- amat[v, ]
			n <- dmat[v, ]

			# Vectorized dbinom over all VAFs x cells in one call:
			# replicate x,n for each VAF (row-major by VAF), replicate p across cells
			X <- rep(x, each = K)
			N <- rep(n, each = K)
			P <- rep(p, times = ncells)

			val <- dbinom(x = X, size = N, prob = P, log = log)

			# Reshape to [K x #cells] so that rows=VAFs, cols=cells
			m <- matrix(val, nrow = K, ncol = ncells, byrow = FALSE)
			rownames(m) <- vafs
			colnames(m) <- cell_names
			m
		}
	)

	names(liks) <- variants
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
get_leaf_liks_mat_old = function(amat, dmat, vafs, eps = 0, ncores = 1, log = FALSE) {

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
                            prob = pmin(vaf + eps, 1 - eps),
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
convert_logliks_to_logP_list_colmajor <- function(logliks, phy) {
	# Build a column-major (C x n) flattened logP_list:
	# for each locus, states (rows) are contiguous per node (column),
	# so the flattened index is node * C + state (0-based).
	E <- reorder_phylo(phy)$edge
	phy$node.label <- NULL

	logP_list_cm <- lapply(logliks, function(liks_mut) {
		n_tips <- length(phy$tip.label)
		n_nodes <- phy$Nnode
		root_node <- E[nrow(E), 1]
		k <- nrow(liks_mut)

		# Matrix shape: C x n (states x nodes); R uses column-major by default.
		P <- matrix(nrow = k, ncol = n_tips + n_nodes)
		rownames(P) <- rownames(liks_mut)

		# Tip log-likelihoods, internal nodes prior, and root node clamp
		P[, 1:n_tips] <- liks_mut[, phy$tip.label]
		P[, (n_tips + 1):(n_tips + n_nodes)] <- log(1 / k)
		P[, root_node] <- c(0, rep(-Inf, k - 1))  # log(c(1, 0, ...))

		# Flatten in column-major order (default in R), yielding node-major layout.
		as.vector(P)
	})

	return(logP_list_cm)
}




#' @export 
get_vaf_bins = function(k) {
    bins = seq(0, 1, 1/k)
    bins = c(0,bins,1)
    vafs = sapply(1:(length(bins)-1), function(i){(bins[i] + bins[i+1])/2})
    return(vafs)
}

#' @export 
get_transition_mat_wf = function(k, eps = 0.01, N = 100, n_rep = 1e4, ngen = 100) {

    A = matrix(NA, ncol = k + 2, k + 2)
    bins = seq(0, 1, 1/k)
    bins = c(0,bins,1)

    for(i in 1:(length(bins)-1)) {

        p = mean(c(bins[i], bins[i+1]))

        p_gen = p
        for (gen in 1:ngen) {
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

    if (ngen == 0) {
        A = diag(k+2)
    }

    vs = sapply(1:(length(bins)-1), function(i){(bins[i] + bins[i+1])/2})
    colnames(A) = vs
    rownames(A) = vs

    return(A)

}


# Caching environments for transition matrices
.mitodrift_T_cache <- new.env(parent = emptyenv())
.mitodrift_Tmat_cache <- new.env(parent = emptyenv())

#' Get the transition matrix for WF model with HMM (with caching)
#' TODO: add log option for small probabilities
#' @param k number of VAF bins
#' @param eps error rate
#' @param N population size
#' @param ngen number of generations
#' @param safe whether to add small probability to avoid 0s
#' @return transition matrix
#' @export
get_transition_mat_wf_hmm <- function(k, eps, N, ngen, safe = FALSE) {
	# Precompute bin boundaries and VAF bin midpoints
	bin_boundaries <- c(0, seq(0, 1, 1 / k), 1)
	vaf_bins <- get_vaf_bins(k)

	# Zero generations: identity (then boundary rows adjusted below)
	if (ngen == 0) {
		A <- diag(k + 2)
	} else {
		# ---- T_ngen cache key depends only on N (matrix size) and ngen (power) ----
		key_T <- paste0(N, "|", as.integer(ngen))

		if (exists(key_T, envir = .mitodrift_T_cache, inherits = FALSE)) {
			T_ngen <- get(key_T, envir = .mitodrift_T_cache, inherits = FALSE)
		} else {
			# Build or reuse the single-generation transition matrix T_mat for this N
			key_Tmat <- as.character(N)
			if (exists(key_Tmat, envir = .mitodrift_Tmat_cache, inherits = FALSE)) {
				T_mat <- get(key_Tmat, envir = .mitodrift_Tmat_cache, inherits = FALSE)
			} else {
				p_vec <- (0:N) / N
				T_mat <- outer(p_vec, 0:N, function(p, k) dbinom(k, size = N, prob = p))
				assign(key_Tmat, T_mat, envir = .mitodrift_Tmat_cache)
			}

			# Multi-generation transition via integer matrix power
			T_ngen <- expm::`%^%`(T_mat, ngen)
			assign(key_T, T_ngen, envir = .mitodrift_T_cache)
		}

		# ---- Aggregate allele-count probabilities into VAF bins ----
		A <- matrix(0, nrow = k + 2, ncol = k + 2)
		for (i in 1:(k + 2)) {
			# Representative starting allele count for bin i (midpoint mapping)
			start_vaf <- vaf_bins[i]
			start_allele_count <- round(start_vaf * N)
			start_allele_count <- max(0, min(N, start_allele_count))

			# Probability distribution after ngen generations
			prob_dist_after_ngen <- T_ngen[start_allele_count + 1, ]

			for (j in 1:(k + 2)) {
				xstart <- floor(bin_boundaries[j] * N) + 1
				xend <- floor(bin_boundaries[j + 1] * N)
				xstart <- min(xstart, xend)
				if (j == k + 1) xend <- xend - 1
				if (xstart <= xend) {
					indices <- xstart:xend
					A[i, j] <- sum(prob_dist_after_ngen[indices + 1])
				} else {
					A[i, j] <- 0
				}
			}
		}
	}

	# ---- Apply mutation rate to boundary rows ----
	A[1, ] <- c(1 - eps, rep(eps / (k + 1), k + 1))
	A[nrow(A), ] <- rev(A[1, ])

	# Names and safety floor
	colnames(A) <- vaf_bins
	rownames(A) <- vaf_bins
	if (safe) {
		A[A == 0] <- 1e-16
	}
	return(A)
}

#' Wrapper function to interpolate non-integer generations
#' @param k number of VAF bins
#' @param eps error rate
#' @param N population size
#' @param ngen number of generations
#' @return transition matrix
#' @export
get_transition_mat_wf_hmm_wrapper = function(k, eps, N, ngen, safe = FALSE) {
    
    # Check if ngen is an integer
    if (ngen == round(ngen)) {
        # If integer, use the original function directly
        return(get_transition_mat_wf_hmm(k = k, eps = eps, N = N, ngen = ngen, safe = safe))
    }
    
    # If non-integer, interpolate between two nearest integer generations
    ngen_lower = floor(ngen)
    ngen_upper = ceiling(ngen)
    
    # Get transition matrices for the two nearest integer generations
    A_lower = get_transition_mat_wf_hmm(k = k, eps = eps, N = N, ngen = ngen_lower, safe = safe)
    A_upper = get_transition_mat_wf_hmm(k = k, eps = eps, N = N, ngen = ngen_upper, safe = safe)
    
    # Calculate interpolation weight
    weight = ngen - ngen_lower
    
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

#' Convert allele count matrix (amat) and total count matrix (dmat) to long format
#' @param amat Allele count matrix with variants as rows and cells as columns
#' @param dmat Total count matrix with variants as rows and cells as columns
#' @return A data.table in long format with columns: variant, cell, a (allele count), d (total count)
#' @export
mat_to_long = function(amat, dmat) {
    # Ensure both matrices have the same dimensions and row/column names
    if (!identical(dim(amat), dim(dmat)) || 
        !identical(rownames(amat), rownames(dmat)) || 
        !identical(colnames(amat), colnames(dmat))) {
        stop("amat and dmat must have identical dimensions and row/column names")
    }
    
    # Convert to data.table and melt to long format
    amat_dt <- as.data.table(amat, keep.rownames = "variant")
    dmat_dt <- as.data.table(dmat, keep.rownames = "variant")
    
    # Melt both matrices to long format
    amat_long <- data.table::melt(amat_dt, id.vars = "variant", 
                                  variable.name = "cell", value.name = "a")
    dmat_long <- data.table::melt(dmat_dt, id.vars = "variant", 
                                  variable.name = "cell", value.name = "d")
    
    # Merge the two long tables
    result <- merge(amat_long, dmat_long, by = c("variant", "cell"))
    
    # Set column order
    setcolorder(result, c("variant", "cell", "a", "d"))
    
    return(result)
}

#' Make a rooted NJ tree
#' @param vmat A matrix of cell-by-variable values
#' @param dist_method The distance method to use
#' @param ncores Number of threads for `parallelDist::parDist` (default: 1)
#' @return A phylo object
#' @export
make_rooted_nj = function(vmat, dist_method = 'manhattan', ncores = 1) {
    vmat[is.na(vmat)] = 0
    vmat = cbind(vmat, outgroup = 0) %>% as.matrix %>% t
    if (ncores > 1) {
        dist_mat = parallelDist::parDist(vmat, method = dist_method, threads = ncores)
    } else {
        dist_mat = dist(vmat, method = dist_method)
    }
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
    res_max = NULL

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
            }
        }

        logZ = c(logZ, res_mar$logZ)

        if (store_bels) {
            ebels[[mut]] = res_mar$edge.bel
            nbels[[mut]] = res_mar$node.bel
            rownames(nbels[[mut]]) = vnames
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
            ngen = flip_df[i,]$length

            A = interpolate_matrices(As, ngen, min_index = 1e-5)

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

    phy_new = phy
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

# Minimal guard: safe read for possibly truncated/partial RDS
safe_read_chain = function(path, ncores = 1) {
    if (!file.exists(path)) return(NULL)
    fi = file.info(path)
    if (is.na(fi$size) || fi$size <= 0) return(NULL)
    tryCatch(qs2::qd_read(path, nthreads = ncores), error = function(e) NULL)
}

# to fix: apparently the init tree has to be rooted otherwise to_phylo_reoder won't work. 
#' @export
run_tree_mcmc_batch = function(
    phy_init, logP_list, logA_vec, outfile, max_iter = 100, nchains = 1, ncores = 1, ncores_qs = 1,
    batch_size = 1000, diag = TRUE, conv_thres = NULL, resume = FALSE
) {


    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    ncores_qs <- if (isTRUE(qs2:::check_TBB())) ncores_qs else 1L
    message('Using ', ncores_qs, ' cores for saving/writing MCMC trace')

    if (!is.null(conv_thres) && !diag) {
        warning('conv_thres provided but diag = FALSE; enabling diagnostics to monitor convergence')
        diag <- TRUE
    }

    chains = 1:nchains

    outdir = dirname(outfile)

    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    edge_list_all = if (resume) safe_read_chain(outfile, ncores = ncores_qs) else NULL
    if (is.null(edge_list_all)) {
        edge_list_all = vector('list', nchains)
    } else {
        length(edge_list_all) = nchains
    }

    max_len = if (is.null(conv_thres)) max_iter + 1L else NULL
    for (i in seq_along(edge_list_all)) {
        chain_list = edge_list_all[[i]]
        if (is.null(chain_list) || length(chain_list) == 0) {
            edge_list_all[[i]] = list(phy_init$edge)
            next
        }
        if (!is.null(max_len) && length(chain_list) > max_len) {
            edge_list_all[[i]] = chain_list[seq_len(max_len)]
        }
    }
    names(edge_list_all) = as.character(chains)

    chain_lengths = vapply(edge_list_all, length, integer(1))
    if (length(unique(chain_lengths)) != 1L) {
        stop('run_tree_mcmc_batch assumes all chains have the same length. Please regenerate or clean the saved state before resuming.')
    }
    completed_iters = chain_lengths[1] - 1L

    if (is.null(conv_thres)) {
        remaining = max_iter - completed_iters
        message('Remaining iterations per chain: ', remaining)
        if (remaining <= 0L) {
            message('All chains have completed the requested iterations.')
            qs2::qd_save(edge_list_all, outfile, nthreads = ncores_qs)
            return(edge_list_all)
        }
        total_batches = ceiling(remaining / batch_size)
        message('Running MCMC with ', length(chains), ' chains in up to ', total_batches, ' batches of ', batch_size)
    } else {
        message('Running MCMC with ', length(chains), ' chains until ASDSF <= ', conv_thres,
                ' (batch size ', batch_size, ')')
    }

    batch_idx = 0L
    repeat {
        chain_lengths = vapply(edge_list_all, length, integer(1))
        if (length(unique(chain_lengths)) != 1L) {
            stop('run_tree_mcmc_batch assumes all chains have the same length. Found inconsistent lengths during execution.')
        }
        completed_iters = chain_lengths[1] - 1L

        if (is.null(conv_thres)) {
            remaining = max_iter - completed_iters
            if (remaining <= 0L) {
                message('All chains have completed the requested iterations.')
                break
            }
            iter_this_batch = min(batch_size, remaining)
            batch_label = paste('batch', batch_idx + 1L, 'of', ceiling(remaining / batch_size))
        } else {
            iter_this_batch = batch_size
            batch_label = paste('batch', batch_idx + 1L)
        }

        batch_idx = batch_idx + 1L
        message('Running ', batch_label)

        ptm <- proc.time()

        iter_vec = rep(iter_this_batch, length(chains))
        start_edges = lapply(edge_list_all, function(chain_list) {
            chain_list[[length(chain_list)]]
        })
        seed_vec = as.integer(1000003L * (batch_idx - 1L) + chains)

        elist_active = tree_mcmc_parallel_seeded(
            start_edges,
            logP_list,
            logA_vec,
            iter_vec,
            seed_vec
        )

        for (chain_id in seq_along(edge_list_all)) {
            elist = restore_elist(elist_active[[chain_id]])
            if (length(elist) > 0) {
                elist = elist[-1]
            }
            edge_list_all[[chain_id]] = c(edge_list_all[[chain_id]], elist)
        }

        qs2::qd_save(edge_list_all, outfile, nthreads = ncores_qs)

        if (diag) {
            asdsf <- compute_target_tree_asdsf(
                phy_target = phy_init,
                edge_list_chains = edge_list_all,
                min_freq = 0,
                rooted = TRUE,
                ncores = ncores
            )
            message('ASDSF (target clades) after ', batch_label, ': ', signif(asdsf, 4))
            if (!is.null(conv_thres) && !is.na(asdsf) && asdsf <= conv_thres) {
                message('Convergence threshold (ASDSF) reached. Stopping MCMC.')
                break
            }
        }

        batch_time <- proc.time() - ptm
        message(paste('Completed', batch_label, paste0('(', signif(batch_time[['elapsed']], 2), 's', ')')))
    }

    qs2::qd_save(edge_list_all, outfile, nthreads = ncores_qs)

    return(edge_list_all)
}

#' Uniform random-walk tree MCMC (accept-all)
#'
#' Mirrors `run_tree_mcmc_batch` but uses the uniform sampler that accepts every NNI proposal
#' (no likelihood evaluation). Useful for generating random walks over tree space quickly.
#' @export
run_tree_mcmc_batch_uniform = function(
    phy_init, outfile, max_iter = 100, nchains = 1, ncores = 1, ncores_qs = 1,
    batch_size = 1000, diag = TRUE, conv_thres = NULL, resume = FALSE
) {


    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    ncores_qs <- if (isTRUE(qs2:::check_TBB())) ncores_qs else 1L
    message('Using ', ncores_qs, ' cores for saving/writing MCMC trace')

    if (!is.null(conv_thres) && !diag) {
        warning('conv_thres provided but diag = FALSE; enabling diagnostics to monitor convergence')
        diag <- TRUE
    }

    chains = 1:nchains

    outdir = dirname(outfile)

    if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }

    edge_list_all = if (resume) safe_read_chain(outfile, ncores = ncores_qs) else NULL
    if (is.null(edge_list_all)) {
        edge_list_all = vector('list', nchains)
    } else {
        length(edge_list_all) = nchains
    }

    max_len = if (is.null(conv_thres)) max_iter + 1L else NULL
    for (i in seq_along(edge_list_all)) {
        chain_list = edge_list_all[[i]]
        if (is.null(chain_list) || length(chain_list) == 0) {
            edge_list_all[[i]] = list(phy_init$edge)
            next
        }
        if (!is.null(max_len) && length(chain_list) > max_len) {
            edge_list_all[[i]] = chain_list[seq_len(max_len)]
        }
    }
    names(edge_list_all) = as.character(chains)

    chain_lengths = vapply(edge_list_all, length, integer(1))
    if (length(unique(chain_lengths)) != 1L) {
        stop('run_tree_mcmc_batch_uniform assumes all chains have the same length. Please regenerate or clean the saved state before resuming.')
    }
    completed_iters = chain_lengths[1] - 1L

    if (is.null(conv_thres)) {
        remaining = max_iter - completed_iters
        message('Remaining iterations per chain: ', remaining)
        if (remaining <= 0L) {
            message('All chains have completed the requested iterations.')
            qs2::qd_save(edge_list_all, outfile, nthreads = ncores_qs)
            return(edge_list_all)
        }
        total_batches = ceiling(remaining / batch_size)
        message('Running uniform MCMC with ', length(chains), ' chains in up to ', total_batches, ' batches of ', batch_size)
    } else {
        message('Running uniform MCMC with ', length(chains), ' chains until ASDSF <= ', conv_thres,
                ' (batch size ', batch_size, ')')
    }

    batch_idx = 0L
    repeat {
        chain_lengths = vapply(edge_list_all, length, integer(1))
        if (length(unique(chain_lengths)) != 1L) {
            stop('run_tree_mcmc_batch_uniform assumes all chains have the same length. Found inconsistent lengths during execution.')
        }
        completed_iters = chain_lengths[1] - 1L

        if (is.null(conv_thres)) {
            remaining = max_iter - completed_iters
            if (remaining <= 0L) {
                message('All chains have completed the requested iterations.')
                break
            }
            iter_this_batch = min(batch_size, remaining)
            batch_label = paste('batch', batch_idx + 1L, 'of', ceiling(remaining / batch_size))
        } else {
            iter_this_batch = batch_size
            batch_label = paste('batch', batch_idx + 1L)
        }

        batch_idx = batch_idx + 1L
        message('Running ', batch_label)

        ptm <- proc.time()

        iter_vec = rep(iter_this_batch, length(chains))
        start_edges = lapply(edge_list_all, function(chain_list) {
            chain_list[[length(chain_list)]]
        })
        seed_vec = as.integer(1000003L * (batch_idx - 1L) + chains)

        elist_active = tree_mcmc_parallel_seeded_uniform(
            start_edges,
            iter_vec,
            seed_vec
        )

        for (chain_id in seq_along(edge_list_all)) {
            elist = restore_elist(elist_active[[chain_id]])
            if (length(elist) > 0) {
                elist = elist[-1]
            }
            edge_list_all[[chain_id]] = c(edge_list_all[[chain_id]], elist)
        }

        qs2::qd_save(edge_list_all, outfile, nthreads = ncores_qs)

        if (diag) {
            asdsf <- compute_target_tree_asdsf(
                phy_target = phy_init,
                edge_list_chains = edge_list_all,
                min_freq = 0,
                rooted = TRUE,
                ncores = ncores
            )
            message('ASDSF (target clades) after ', batch_label, ': ', signif(asdsf, 4))
            if (!is.null(conv_thres) && !is.na(asdsf) && asdsf <= conv_thres) {
                message('Convergence threshold (ASDSF) reached. Stopping uniform MCMC.')
                break
            }
        }

        batch_time <- proc.time() - ptm
        message(paste('Completed', batch_label, paste0('(', signif(batch_time[['elapsed']], 2), 's', ')')))
    }

    qs2::qd_save(edge_list_all, outfile, nthreads = ncores_qs)

    return(edge_list_all)
}

restore_elist = function(elist) {
    lapply(elist, function(edges){matrix(edges, ncol = 2)})
}

#' @export
collect_chains = function(edge_list_all, phy_init, burnin = 0, max_iter = Inf) {

    if (max_iter < burnin) {
        stop('Max iter needs to be greater than burnin')
    }

    # drop empty chains
    edge_list_all = edge_list_all[edge_list_all %>% sapply(length) > 0]

    # reconstruct trees per chain and apply burnin/truncation
    mcmc_trees = edge_list_all %>% lapply(function(elist){
            elist = elist[(burnin+1):min(length(elist), max_iter)]
            lapply(elist, function(edges){
                tree = attach_edges(phy_init, edges)
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

collect_edges = function(edge_list_all, burnin = 0, max_iter = Inf) {

    if (max_iter < burnin) {
        stop('Max iter needs to be greater than burnin')
    }

    mcmc_edges = edge_list_all %>% lapply(function(elist){
            elist = elist[(burnin+1):min(length(elist), max_iter)]
        }) %>% unlist(recursive = F)

    return(mcmc_edges)
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


#' Parallelized Clade Support Calculation
#'
#' Computes the support for each clade in a reference phylogeny (`phy`) across a list of phylogenetic trees (`phy_list`).
#' The calculation can be parallelized across multiple cores for efficiency.
#'
#' @param phy A reference phylogeny of class \code{phylo}.
#' @param phy_list A list of phylogenetic trees (of class \code{phylo}) to compare against the reference tree.
#' @param rooted Logical; whether to treat the trees as rooted. Default is \code{FALSE}.
#' @param ncores Integer; number of cores to use for parallel computation. Default is \code{1} (no parallelization).
#' @param normalize Logical; if \code{TRUE}, the support values are normalized to the range [0, 1] by dividing by the number of trees. Default is \code{TRUE}.
#'
#' @return A numeric vector of clade support values, in the same order as the internal nodes of \code{phy}.
#'         If \code{normalize = TRUE}, values are proportions; otherwise, they are counts.
#'
#' @export
prop.clades.par <- function(phy, phy_list, rooted = FALSE,
                            ncores = 1, 
                            normalize = TRUE) {

    nT <- length(phy_list)
    if (nT == 0) stop("Need at least one tree in the tree list")

    # If ncores = 1, skip chunking and parallel processing
    if (ncores == 1) {
        support <- prop.clades(phy, phy_list, rooted = rooted)
    } else {
        # split indices into mc.cores chunks
        chunks <- split(seq_len(nT), cut(seq_len(nT), breaks = ncores, labels = FALSE))

        # for each chunk, compute support counts
        partials <- mclapply(chunks, function(idxs) {
            pc = prop.clades(phy, phy_list[idxs], rooted = rooted)
            pc[is.na(pc)] = 0
            return(pc)
        }, mc.cores = ncores)

        # sum support across chunks
        support <- Reduce(`+`, partials)
    }

    if (normalize) {
        support <- support / nT
    }

    support[is.na(support)] = 0

    return(support)
}

#' Add clade frequencies to a phylogenetic tree
#'
#' Computes the frequency (support) of each clade in a reference phylogeny (`phy`) across a list of phylogenetic trees (`edge_list`),
#' and adds these frequencies as node labels to the reference tree.
#'
#' @param phy A reference phylogeny of class \code{phylo}.
#' @param edge_list A list of edge matrices (or phylogenetic trees) to compare against the reference tree.
#' @param rooted Logical; whether to treat the trees as rooted. Default is \code{TRUE}.
#' @param ncores Integer; number of cores to use for parallel computation. Default is \code{1} (no parallelization).
#' @return A phylo object with clade frequencies added as node labels.
#' @export
add_clade_freq = function(phy, edge_list, rooted = TRUE, ncores = 1) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)
    phy = reorder_phylo(phy)
    freqs = prop_clades_par(phy$edge, edge_list, rooted = rooted, normalize = TRUE)
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
    if (!is.binary(tree)) {stop("Tree must be binary")}
    node_confs = as.numeric(tree$node.label)
    node_confs[is.na(node_confs)] = 0
    tree = TreeTools::CollapseNode(tree, which(node_confs < conf) + length(tree$tip.label))
    tree = TreeTools::Renumber(tree)
    return(tree)
}


#' @export
trim_tree_size = function(tree, min_conf = 0, min_frac = 0, max_frac = Inf, method = 'union') {
    node_confs = as.numeric(tree$node.label)
    node_confs[is.na(node_confs)] = 0
    
    ntip = length(tree$tip.label)
    node_ids = (ntip + 1):(ntip + tree$Nnode)
    max_size = ceiling(ntip * max_frac)
    min_size = ceiling(ntip * min_frac)
    
    # Calculate clade sizes
    if (requireNamespace("phangorn", quietly = TRUE)) {
         sizes = lengths(phangorn::Descendants(tree, node_ids, type = "tips"))
    } else {
         sizes = sapply(node_ids, function(x) length(ape::extract.clade(tree, x)$tip.label))
    }
    
    if (method == 'union') {
        # Collapse if EITHER criteria is failed (strict filtering)
        # i.e. keep only if conf >= min_conf AND size >= min_size AND size <= max_size
        to_collapse = which(node_confs < min_conf | sizes < min_size | sizes > max_size)
    } else if (method == 'intersection') {
        # Collapse only if BOTH criteria are failed (lenient filtering)
        # i.e. keep if conf >= min_conf OR (size >= min_size AND size <= max_size)
        to_collapse = which(node_confs < min_conf & (sizes < min_size | sizes > max_size))
    } else {
        stop("method must be 'union' or 'intersection'")
    }
    
    if (length(to_collapse) > 0) {
        tree = TreeTools::CollapseNode(tree, node_ids[to_collapse])
        tree = TreeTools::Renumber(tree)
    }
    return(tree)
}

#' @export
trim_tree_exp = function(tree, tol) {
	
	n_tip <- length(tree$tip.label)
    n_tol = n_tip * tol
	
	# branch posteriors on internal nodes
	node_confs <- as.numeric(tree$node.label)
	node_confs[is.na(node_confs)] <- 0
	
	internal_nodes <- (n_tip + 1L):(n_tip + tree$Nnode)
	desc_nodes <- TreeTools::CladeSizes(tree, nodes = internal_nodes)
	clade_tips <- as.integer((desc_nodes + 2L) / 2L)
	
	# Expected wrong cells per branch: (1 - p) * clade_size
	exp_wrong <- (1 - node_confs) * clade_tips

	# Collapse branches whose expected wrong cells exceeds epsilon (= n_tol)
	nodes_to_collapse <- which(exp_wrong > n_tol) + n_tip
	
	if (length(nodes_to_collapse) > 0L) {
		tree <- TreeTools::CollapseNode(tree, nodes_to_collapse)
		tree <- TreeTools::Renumber(tree)
	}
	
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

#' Add "Node<n>" labels to the internal nodes of a phylo tree
#'
#' @param tree A phylo object
#' @param prefix Character prefix for node names (default "Node")
#' @return The same phylo object, with tree$node.label set to prefix + node numbers
add_node_names <- function(tree, prefix = "Node", start_from_tip = TRUE) {
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylo object")
  }
  ntip  <- length(tree$tip.label)
  nnode <- tree$Nnode

  # internal node IDs run from ntip+1 to ntip+nnode
  if (start_from_tip) {
    ids <- seq(ntip + 1, ntip + nnode)
  } else {
    ids <- seq(1, nnode)
  }

  # assign labels
  tree$node.label <- paste0(prefix, ids)

  return(tree)
}

map_cell_to_tree = function(tree, cell, logliks, logA_vec, leaf_only = FALSE, ncores = 1) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    tree = TreeTools::Renumber(tree)
    
    if (leaf_only) {
        nodes = 1:length(tree$tip.label)
        labels = tree$tip.label
    } else {
        nodes = 1:(length(tree$tip.label) + tree$Nnode)
        tree = add_node_names(tree)
        labels = c(tree$tip.label, tree$node.label)
    }

    # check if all cells are in logliks
    if (!all(c(cell, tree$tip.label) %in% colnames(logliks[[1]]))) {
        stop(glue('Some cells not found in logliks'))
    }

    # generate all possible assignment to nodes
    trees_new = lapply(
        nodes,
        function(i) {
            TreeTools::AddTip(tree, where = i, label = cell)
        })

    edges_new = lapply(
        trees_new, 
        function(tree) {
            reorderRcpp(tree$edge)
        })

    # subset likelihoods
    logP_list = convert_logliks_to_logP_list(logliks, trees_new[[1]])

    # score assignments
    scores = score_trees_parallel(edges_new, logP_list, logA_vec)
    
    probs = exp(scores - logSumExp(scores))
    
    res = data.frame(guide_node = labels, cell_map = cell, probs = probs, scores = scores)

    return(res)
}

save_qd_par = function(objects, paths, ncores = 1) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    save_qd_cpp(objects, paths)
}