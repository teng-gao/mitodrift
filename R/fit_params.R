#' Fit tree parameters using EM (BP-backed, native ordering)
#'
#' A faster EM that uses compute_node_edge_beliefs_bp2() directly.
#' It relies on the inherent ordering used by the C++ routine:
#'   - node_beliefs[[locus]] is an n x C matrix in APE node-id order (1..n)
#'   - edge_beliefs[[locus]] is an m x C x C array where the m edges correspond
#'     to the post-ordered edge list internal to the C++ routine (order is irrelevant
#'     because we sum over edges).
#'
#' @param tree_fit phylogenetic tree (will be renumbered)
#' @param amat alternative allele counts
#' @param dmat total depth matrix
#' @param initial_params named numeric: ngen, log_eps, log_err
#' @param lower_bounds named numeric lower bounds (same names)
#' @param upper_bounds named numeric upper bounds (same names)
#' @param max_iter maximum EM iterations
#' @param k number of hidden states (VAF bins)
#' @param npop population size for WF-HMM transition
#' @param ncores kept for API compatibility (unused here)
#' @param epsilon convergence threshold on parameter deltas
#' @param trace if TRUE, returns trace dataframe
#' @return either named vector of final params (ngen, eps, err) or list(par=..., trace=...)
#' @keywords internal
fit_params_em_cpp <- function(tree_fit, amat, dmat,
	initial_params = c('ngen' = 100, 'log_eps' = log(1e-3), 'log_err' = log(1e-3)),
	lower_bounds = c('ngen' = 1, 'log_eps' = log(1e-12), 'log_err' = log(1e-12)),
	upper_bounds = c('ngen' = 1000, 'log_eps' = log(0.2), 'log_err' = log(0.2)),
	max_iter = 10, k = 20, npop = 600, ncores = 1, epsilon = 1e-3, trace = TRUE) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

	# Ensure stable numeric node IDs 1..n (matches C++ expectations)
	tree_fit <- TreeTools::Renumber(tree_fit)
	par <- initial_params
	if (trace) trace_df <- data.frame()

	# Precompute VAF grid
	vafs <- get_vaf_bins(k)

	for (i in seq_len(max_iter)) {
		ngen <- par[['ngen']]
		eps  <- exp(par[['log_eps']])
		err  <- exp(par[['log_err']])

		# ---------- E step (BP, native ordering) ----------
		# Emissions in log-space, list of (C x |S_l|) matrices with colnames as node ids
		logliks <- get_leaf_liks_mat_cpp(amat, dmat, vafs = vafs, eps = err, log = TRUE)
		# Pack to C*n vectors per locus in APE 1..n order
		logP_list <- convert_logliks_to_logP_list(logliks, tree_fit)

		# Transition matrix (row=parent state, col=child state), row-major packing
		A <- get_transition_mat_wf_hmm(k = k, eps = eps, N = npop, ngen = ngen, safe = TRUE)
		logA_vec <- as.vector(t(log(A)))

		# Call the streamlined BP stats routine
		res_cpp <- compute_node_edge_stats_bp2(tree_fit$edge, logP_list, logA_vec)
		node_post <- res_cpp$leaf_beliefs
		E_counts <- res_cpp$edge_counts
		logL <- res_cpp$logZ

		ids <- match(colnames(logliks[[1]]), tree_fit$tip.label)
		leaf_post <- lapply(node_post, function(x) x[ids, , drop = FALSE])

		message(glue("Iteration {i} E step: logL = {logL}"))
		if (trace) {
			trace_df <- dplyr::bind_rows(trace_df, data.frame(
				iter = i,
				ngen = par[['ngen']],
				log_eps = par[['log_eps']],
				log_err = par[['log_err']],
				logL = logL
			))
		}

		# ---------- M step (split) ----------
		# 1) Optimize (ngen, log_eps) using ONLY the transition term
		Q_trans <- function(x) {
			A_cur <- get_transition_mat_wf_hmm_wrapper(k = k, eps = exp(x[2]), N = npop, ngen = x[1], safe = TRUE)
			-sum(E_counts * log(A_cur))
		}

		fit_trans <- optim(
			par = c(par[['ngen']], par[['log_eps']]),
			fn = Q_trans, method = 'L-BFGS-B',
			lower = c(lower_bounds[['ngen']], lower_bounds[['log_eps']]),
			upper = c(upper_bounds[['ngen']], upper_bounds[['log_eps']])
		)

		# 2) Optimize log_err using ONLY the emission term
		Q_leaf <- function(le) {
			err_cur <- exp(le)
			logliks <- get_leaf_liks_mat_cpp(amat, dmat, vafs = vafs, eps = err_cur, log = TRUE)
			sum(vapply(seq_along(logliks), function(ii) {
				sum(leaf_post[[ii]] * t(logliks[[ii]]))
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
		if (i > 1 && all(abs(par - par_new) < epsilon)) {
			message(glue("Converged at iteration {i}"))
			par <- par_new
			break
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