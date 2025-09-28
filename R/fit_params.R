#' Fit tree parameters using EM algorithm
#' @param tree_fit phylogenetic tree
#' @param amat alternative allele count matrix
#' @param dmat total depth matrix
#' @param initial_params initial parameter values (ngen, log_eps, log_err)
#' @param lower_bounds lower bounds for parameters in transformed space
#' @param upper_bounds upper bounds for parameters in transformed space
#' @param max_iter maximum number of iterations
#' @param k number of clusters for likelihood computation
#' @param npop population size for likelihood computation
#' @param ncores number of cores to use
#' @param epsilon convergence threshold
#' @param trace whether to return trace of parameter values
#' @return parameter values or list of parameter values and trace
#' @export
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


#' Decode tree using EM algorithm
#' @param tn phylogenetic tree
#' @param A transition matrix
#' @param liks leaf likelihoods
#' @param ncores number of cores to use
#' @return list of logZ, edge beliefs, and node beliefs
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



#' Compute the likelihood of a tree given model parameters
#' @param tree_fit phylogenetic tree
#' @param ngen number of generations
#' @param err error rate
#' @param eps mutation rate
#' @param npop population size
#' @return log-likelihood
#' @export
get_param_lik_cpp = function(tree_fit, amat, dmat, ngen, err, eps, npop = 600, k = 20) {
    
    A = get_transition_mat_wf_hmm(k = k, eps = eps, N = npop, ngen = ngen)
    logA_vec = t(log(A))
    
    logliks = get_leaf_liks_mat(amat, dmat, get_vaf_bins(k = k), eps = err, log = TRUE)
    logP_list = convert_logliks_to_logP_list(logliks, tree_fit)

    l = mitodrift:::score_tree_bp_wrapper2(tree_fit$edge, logP_list = logP_list, logA = logA_vec)

    return(l)
}



#' Fit tree parameters using MCMC
#' 
#' Estimates tree parameters (number of generations, error rates) using MCMC sampling
#' with the fmcmc package.
#' 
#' @param tree_fit phylogenetic tree
#' @param amat alternative allele count matrix
#' @param dmat total depth matrix
#' @param initial_params initial parameter values (ngen, log_eps, log_err)
#' @param lower_bounds lower bounds for parameters in transformed space
#' @param upper_bounds upper bounds for parameters in transformed space
#' @param nsteps number of MCMC steps
#' @param nchains number of MCMC chains
#' @param npop population size for likelihood computation
#' @param k number of clusters for likelihood computation
#' @param ncores number of cores to use
#' @param keep number of MCMC steps to keep at the end of each chain
#' @param outfile file to save MCMC result
#' @param check_conv whether to check convergence of parameter fitting
#' @return MCMC result object from fmcmc
#' @export
fit_params_mcmc = function(
    tree_fit, amat, dmat,
    initial_params = c('ngen' = 100, 'log_eps' = log(1e-3), 'log_err' = log(1e-3)),
    lower_bounds = c('ngen' = 1, 'log_eps' = log(1e-18), 'log_err' = log(1e-18)),
    upper_bounds = c('ngen' = 1000, 'log_eps' = log(0.2), 'log_err' = log(0.2)),
    nsteps = 500, nchains = 1, outfile = NULL,
    npop = 600, k = 20, ncores = 1, keep = 100, check_conv = FALSE) {
    
    # Ensure tree is properly formatted
    tree_fit = reorder_phylo(tree_fit)
    tree_fit$edge.length = NULL
    
    # All three parameters will be estimated
    param_names = c('ngen', 'log_eps', 'log_err')
    
    message(glue("Estimating all 3 parameters using MCMC: {paste(param_names, collapse=', ')}"))
    
    # Create log-likelihood function for MCMC
    # This function takes parameters in transformed space and converts them back
    log_likelihood = function(p, tree_fit., amat., dmat., npop., k.) {
        # Convert parameters back to original scale
        params_orig = p
        
        # Convert log-scale parameters back to original scale
        params_orig['log_eps'] = exp(p['log_eps'])
        names(params_orig)[names(params_orig) == 'log_eps'] = 'eps'
        params_orig['log_err'] = exp(p['log_err'])
        names(params_orig)[names(params_orig) == 'log_err'] = 'err'
        
        # Compute likelihood
        ll = get_param_lik_cpp(tree_fit., amat., dmat., 
            ngen = params_orig['ngen'],
            eps = params_orig['eps'], 
            err = params_orig['err'],
            npop = npop., k = k.)
        
        return(ll)
    }
        
    # Run MCMC
    ncores = min(ncores, nchains)
    message("Starting MCMC sampling using ", ncores, " cores")

    cl <- makeCluster(ncores)

    msg = clusterEvalQ(cl, {
        library(mitodrift)
        library(parallel)
    })

    if (check_conv) {
        checker = fmcmc::convergence_gelman(freq = 50, threshold = 1)
    } else {
        checker = NULL
    }

    kernel = fmcmc::kernel_normal_reflective(
        scale = c(5, 0.1, 0.1),
        scheme = "joint",
        lb = lower_bounds,
        ub = upper_bounds
    )

    mcmc_result = fmcmc::MCMC(
        log_likelihood,
        initial = initial_params,
        nsteps = nsteps,
        nchains = nchains,
        kernel = kernel,
        tree_fit. = tree_fit,
        amat. = amat,
        dmat. = dmat,
        npop. = npop,
        k. = k,
        cl = cl,
        multicore = TRUE,
        conv_checker = checker
    )

    stopCluster(cl)

    res_df_all = lapply(
        seq_len(nchains),
        function(i) {
            data.frame(mcmc_result[[i]]) %>% 
                tibble::rowid_to_column('iter') %>%
                mutate(chain = i)
        }) %>%
        bind_rows() %>%
        mutate(
            err = exp(log_err),
            eps = exp(log_eps),
        ) %>%
        reshape2::melt(id.var = c('iter', 'chain'))

    if (!is.null(outfile)) {
        fwrite(res_df_all, outfile)
    }

    return(res_df_all)
}
