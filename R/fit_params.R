fit_params_em = function(tree_fit, amat, dmat, 
	initial_params = c('ngen' = 100, 'log_eps' = log(1e-3), 'log_err' = log(1e-3)),
	lower_bounds = c('ngen' = 1, 'log_eps' = log(1e-12), 'log_err' = log(1e-12)),
	upper_bounds = c('ngen' = 1000, 'log_eps' = log(0.2), 'log_err' = log(0.2)),
	max_iter = 10, k = 20, npop = 600, ncores = 1, epsilon = 1e-3, trace = TRUE
) {

	# always renumber tree otherwise results are different
	tree_fit <- TreeTools::Renumber(tree_fit)

	par <- initial_params

	if (trace) {
		trace_df <- data.frame()
	}

	# precompute VAF grid once
	vafs <- get_vaf_bins(k)

	for (i in 1:max_iter) {

		ngen <- par[['ngen']]
		eps  <- exp(par[['log_eps']])
		err  <- exp(par[['log_err']])

		# ---------- E step ----------
		liks <- get_leaf_liks_mat_cpp(amat, dmat, vafs = vafs, eps = err)
		A <- get_transition_mat_wf_hmm(k = k, eps = eps, N = npop, ngen = ngen, safe = TRUE)

		res <- decode_tree_em(tree_fit, A, liks, ncores = ncores)
		nbels <- res$nbels
		ebels <- res$ebels
		logL <- res$logZ

		message(glue("Iteration {i} E step: logL = {logL}"))

		if (trace) {
			trace_df <- bind_rows(trace_df, data.frame(
				iter = i,
				ngen = par[['ngen']],
				log_eps = par[['log_eps']],
				log_err = par[['log_err']], 
				logL = logL
			))
		}

		# Sufficient stats for transitions: sum over edges and mutations
		# E_counts is C x C matrix of expected transition counts
		E_counts <- Reduce(`+`, lapply(ebels, function(el) Reduce(`+`, el)))

		# ---------- M step (split) ----------
		# 1) Optimize (ngen, log_eps) using ONLY the transition term:
		#    Q_trans(ngen, log_eps) = sum(E_counts * log A(ngen, eps))
		Q_trans <- function(x) {
			A_cur <- get_transition_mat_wf_hmm_wrapper(k = k, eps = exp(x[2]), N = npop, ngen = x[1], safe = TRUE)
			-sum(E_counts * log(A_cur))
		}

		fit_trans <- optim(
			par = c(par[['ngen']], par[['log_eps']]),
			fn = Q_trans,
			method = 'L-BFGS-B',
			lower = c(lower_bounds[['ngen']], lower_bounds[['log_eps']]),
			upper = c(upper_bounds[['ngen']], upper_bounds[['log_eps']])
		)

		# 2) Optimize log_err using ONLY the emission term:
		#    Q_leaf(log_err) = sum_i sum_{cells,states} nbels_i * logliks_i(err)

		Q_leaf <- function(le) {
			err_cur <- exp(le)
			logliks <- get_leaf_liks_mat_cpp(amat, dmat, vafs = vafs, eps = err_cur, log = TRUE)
			sum(vapply(seq_along(nbels), function(ii) {
				leafs <- colnames(logliks[[ii]])
				sum(nbels[[ii]][leafs, ] * t(logliks[[ii]]))
			}, numeric(1)))
		}

		fit_err <- optimize(
			f = function(le) -Q_leaf(le),
			interval = c(lower_bounds[['log_err']], upper_bounds[['log_err']]),
			maximum = FALSE
		)

		par_new <- c(
			'ngen' = fit_trans$par[1],
			'log_eps' = fit_trans$par[2],
			'log_err' = fit_err$minimum
		)

		# Convergence check
		if (i > 1) {
			if (all(abs(par - par_new) < epsilon)) {
				message(glue("Converged at iteration {i}"))
				par <- par_new
				break
			}
		}

		par <- par_new
		message(glue("Iteration {i} M step: ngen = {par[1]}, log_eps = {par[2]}, log_err = {par[3]}"))
	}

	par_final <- c('ngen' = par[['ngen']], 'eps' = exp(par[['log_eps']]), 'err' = exp(par[['log_err']]))

	if (trace) {
		return(list(par = par_final, trace = trace_df))
	} else {
		return(par_final)
	}
}


#' @export 
decode_tree_em = function(
    tn, A, liks, ncores = 1
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

    res_list <- parallel::mclapply(
        names(liks),
        mc.cores = ncores,
        function(mut) {
            crf_local <- copy_crf(crf)
            # fill missing nodes with uniform prior
            crf_local$node.pot[,] <- 1 / k
            # plug in observed node likelihoods for this mutation
            nodes <- colnames(liks[[mut]])
            if (length(nodes)) crf_local$node.pot[nodes, ] <- t(liks[[mut]])
            # enforce root prior
            crf_local$node.pot[root_node, ] <- c(1, rep(0, k - 1))
            # run exact tree inference
            res_mar <- infer.tree(crf_local)
            rownames(res_mar$node.bel) <- vnames
            list(
                mut = mut,
                logZ = res_mar$logZ,
                edge.bel = res_mar$edge.bel,
                node.bel = res_mar$node.bel
            )
        }
    )

    # assemble outputs exactly like before
    logZ  <- setNames(vapply(res_list, `[[`, numeric(1), "logZ"),
                    vapply(res_list, `[[`, character(1), "mut"))
    ebels <- setNames(lapply(res_list, `[[`, "edge.bel"), names(logZ))
    nbels <- setNames(lapply(res_list, `[[`, "node.bel"), names(logZ))

    logZ = setNames(logZ, names(liks)) 

    return(list('logZ' = sum(logZ), 'ebels' = ebels, 'nbels' = nbels))
}
