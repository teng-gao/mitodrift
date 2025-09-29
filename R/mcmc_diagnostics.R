## -----------------------------------------------------------------------------
#' Compute ASDSF across chains using clades from a target tree
#'
#' @param phy_target A rooted `phylo` object defining the reference clades.
#' @param edge_list_chains List of chains, each a list of edge matrices (2-column integer matrices).
#' @param rooted Logical; treat trees as rooted when matching clades. Default `TRUE`.
#' @param normalize Logical; pass through to `prop_clades_par` (default `TRUE`).
#' @param min_freq Minimum clade frequency threshold used when averaging SDs. Default `0.1`.
#'
#' @return A list with elements `asdsf`, `per_clade_sd`, `keep_mask`, and `freq_matrix`.
#' @export
compute_target_tree_asdsf <- function(phy_target,
                                      edge_list_chains,
                                      rooted = TRUE,
                                      normalize = TRUE,
									  ncores = 1,
                                      min_freq = 0) {

    if (!inherits(phy_target, "phylo")) {
        stop('`phy_target` must be a phylo object')
    }
    if (!is.list(edge_list_chains) || length(edge_list_chains) == 0) {
        stop('`edge_list_chains` must be a non-empty list of chains')
    }

	RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    RcppParallel::setThreadOptions(numThreads = ncores)

    # Compute clade frequencies per chain relative to the target tree
    freqs <- lapply(edge_list_chains, function(chain_edges) {
        if (!length(chain_edges)) {
            # No samples yet: return zeros for all target clades
            return(rep(0, phy_target$Nnode))
        }
        prop_clades_par(
            E_target = phy_target$edge,
            edges = chain_edges,
            rooted = rooted,
            normalize = normalize
        )
    })

    # Ensure all chains yield vectors of identical length
    K <- length(freqs[[1]])
    freq_matrix <- do.call(rbind, freqs)

    if (ncol(freq_matrix) != K) {
        stop('Mismatch in clade count across chains when computing ASDSF')
    }

    # Standard deviation per clade across chains (handle single-chain case)
    if (nrow(freq_matrix) <= 1L) {
        per_clade_sd <- rep(0, K)
    } else {
        per_clade_sd <- apply(freq_matrix, 2, stats::sd)
        per_clade_sd[is.na(per_clade_sd)] <- 0
    }

    keep_mask <- apply(freq_matrix, 2, function(vals) any(vals > min_freq))
    if (!any(keep_mask)) {
        asdsf <- 0
    } else {
        asdsf <- mean(per_clade_sd[keep_mask])
    }

    return(asdsf)
}

#' Sliding-window ASDSF for a target tree
#'
#' Applies `compute_target_tree_asdsf()` to sequential, fixed-width windows of
#' each MCMC chain and returns an ASDSF trajectory where each point uses all
#' samples up to the end of the corresponding window (cumulative windows, akin to
#' `get_asdsfs`).
#'
#' @param phy_target A rooted `phylo` object defining the reference clades.
#' @param edge_list_chains List of chains, each a list of edge matrices (2-column
#'   integer matrices) ordered by iteration.
#' @param window_size Positive integer; number of successive samples in the base
#'   window. Cumulative estimates are computed at multiples of this size.
#' @param burnin Non-negative integer; number of initial samples to discard from
#'   each chain before windowing.
#' @param rooted Logical; treat trees as rooted when matching clades (passed to
#'   `compute_target_tree_asdsf`).
#' @param normalize Logical; passed to `compute_target_tree_asdsf`.
#' @param min_freq Minimum clade frequency threshold when averaging SDs. Passed
#'   to `compute_target_tree_asdsf`.
#' @param min_chains Minimum number of chains required after filtering.
#'
#' @return A data.frame with one row per window containing the ASDSF value, the
#'   window bounds (1-indexed, post-burnin), and attributes `window_size`, `step`,
#'   and `burnin`.
#'
#' @seealso [compute_target_tree_asdsf]
#' @export
compute_target_tree_asdsf_sliding <- function(phy_target,
                                             edge_list_chains,
                                             window_size = 50L,
                                             burnin = 0L,
                                             rooted = TRUE,
                                             normalize = TRUE,
                                             min_freq = 0,
											 ncores = 1,
                                             min_chains = 1L) {

    window_size <- as.integer(window_size)
    burnin <- as.integer(burnin)
    min_chains <- as.integer(min_chains)

    # Drop burn-in samples from each chain
    chains_trim <- lapply(edge_list_chains, function(chain) {
        if (!length(chain)) return(list())
        n <- length(chain)
        if (burnin >= n) {
            list()
        } else {
            chain[(burnin + 1L):n]
        }
    })

    len_vec <- vapply(chains_trim, length, integer(1))
    keep_idx <- which(len_vec >= window_size)
    if (!length(keep_idx)) {
        warning('No chains have at least `window_size` samples after burn-in; returning empty result.')
        return(data.frame())
    }
    if (length(keep_idx) < length(chains_trim)) {
        chains_trim <- chains_trim[keep_idx]
        len_vec <- len_vec[keep_idx]
    }

    min_len <- min(len_vec)
    if (min_len < window_size) {
        warning('Minimum chain length is smaller than `window_size`; returning empty result.')
        return(data.frame())
    }

    ends <- seq.int(window_size, min_len, by = window_size)
    if (!length(ends)) {
        warning('No windows available with the requested parameters; returning empty result.')
        return(data.frame())
    }

    starts <- pmax(1L, ends - window_size + 1L)

    asdsf_vals <- numeric(length(ends))
    for (idx in seq_along(ends)) {
        end <- ends[idx]
        window_chains <- lapply(chains_trim, function(chain) chain[seq_len(end)])
        asdsf_vals[idx] <- compute_target_tree_asdsf(
            phy_target = phy_target,
            edge_list_chains = window_chains,
            rooted = rooted,
            normalize = normalize,
			ncores = ncores,
            min_freq = min_freq
        )
    }

    summary_df <- data.frame(
        window = seq_along(ends),
        window_start = starts + burnin,
        window_end = ends + burnin,
        cumulative_samples = ends,
        n_chains = length(chains_trim),
        window_size = window_size,
        asdsf = asdsf_vals,
        stringsAsFactors = FALSE
    )

    attr(summary_df, 'window_size') <- window_size
    attr(summary_df, 'burnin') <- burnin

    summary_df
}


#' Optimized version of `get_slide_freq_table` from rwty
get_slide_freq_table = function(tree.list, burnin = 0, window.size = 20, gens.per.tree = 1) {
	# Optimized: avoid per-window merges; build a dense matrix of clade frequencies.
	# Keeps original behavior: uses only FULL windows after burnin; if none, returns empty data.frame.
	if (burnin < 0L) burnin <- 0L
	tree.list <- tree.list[(burnin + 1L):length(tree.list)]
	N <- length(tree.list)
	if (!N) return(data.frame())
	n.windows <- as.integer(N / window.size)
	if (n.windows < 1L) return(data.frame())
	# slice into full windows
	idx <- seq_len(n.windows * window.size)
	tree.windows <- split(tree.list[idx], ceiling(seq_along(idx) / window.size))
	# compute clade frequencies per window using existing clade.freq()
	clade.freq.list <- lapply(tree.windows, rwty::clade.freq, start = 1, end = window.size)
	win_names <- prettyNum(seq((burnin + 1L), (burnin + length(clade.freq.list) * window.size), by = window.size) * gens.per.tree, sci = TRUE)
	names(clade.freq.list) <- win_names
	# union of clades across windows; fill a dense matrix (rows = clades, cols = windows)
	all_keys <- unique(unlist(lapply(clade.freq.list, function(df) df$cladenames), use.names = FALSE))
	if (!length(all_keys)) return(data.frame())
	mat <- matrix(0, nrow = length(all_keys), ncol = length(clade.freq.list), dimnames = list(all_keys, win_names))
	for (w in seq_along(clade.freq.list)){
		df <- clade.freq.list[[w]]
		if (!nrow(df)) next
		mat[df$cladenames, w] <- df$cladefreqs
	}
	# convert to data.frame and attach short ids + key map
	df <- as.data.frame(mat, stringsAsFactors = FALSE, optional = TRUE)
	short_ids <- shorten_clade_ids(rownames(df))
	clade_map <- data.frame(id = short_ids, key = rownames(df), stringsAsFactors = FALSE)
	rownames(df) <- short_ids
	attr(df, "clade_key_map") <- clade_map
	# summary columns and ordering (same as original)
	# use matrixStats for speed
	df$sd <- matrixStats::rowSds(as.matrix(df[, win_names, drop = FALSE]))
	df$mean <- matrixStats::rowMeans2(as.matrix(df[, win_names, drop = FALSE]))
	df <- df[order(df$sd, decreasing = TRUE), , drop = FALSE]
	df
}

#' Optimized version of `get_asdsfs` from rwty
get_asdsfs = function(slide.freq.list, min.freq = 0.1){
	# Fast ASDSF over cumulative windows.
	# slide.freq.list: list of slide.freq.tables (data.frames) with window columns + trailing 'sd','mean'.
	stopifnot(is.list(slide.freq.list), length(slide.freq.list) >= 2L)
	x <- slide.freq.list
	sets <- length(x)

	# Identify window columns once
	col_exclude <- intersect(c("sd", "mean"), colnames(x[[1]]))
	win_names <- setdiff(colnames(x[[1]]), col_exclude)
	if (!length(win_names)) stop("No window columns found in slide.freq.table")
	W <- length(win_names)

	# Global clade set (use existing short rownames)
	clades <- unique(unlist(lapply(x, rownames), use.names = FALSE))
	K <- length(clades)

	# Use matrixStats functions directly
	row_cumsums <- matrixStats::rowCumsums
	row_sds <- matrixStats::rowSds
	row_anys <- matrixStats::rowAnys

	# Precompute per-chain cumulative sums (K_j x W) and row index maps into global clades
	CS_list <- vector("list", sets)
	IDX_list <- vector("list", sets)
	for (j in seq_len(sets)){
		# window matrix for chain j
		Mj <- as.matrix(x[[j]][, win_names, drop = FALSE])
		# map rows into global clade order
		idx <- match(rownames(x[[j]]), clades)
		IDX_list[[j]] <- idx
		# cumulative sums across windows per clade (rows)
		CSj <- row_cumsums(Mj)
		CS_list[[j]] <- CSj
	}

	# Allocate result vector
	ASDSF <- numeric(W)

	# Work matrix (K x sets), reused across windows
	A <- matrix(0, nrow = K, ncol = sets)

	for (i in seq_len(W)){
		# fill A with cumulative means at window i
		# zero once per iteration
		A[] <- 0
		inv_i <- 1.0 / i
		for (j in seq_len(sets)){
			v <- CS_list[[j]][, i] * inv_i
			A[IDX_list[[j]], j] <- v
		}
		# min.freq filter at this step
		keep <- row_anys(A > min.freq)
		if (!any(keep)){
			ASDSF[i] <- 0
			next
		}
		per_sd <- row_sds(A[keep, , drop = FALSE])
		ASDSF[i] <- mean(per_sd)
	}

	# Generations from window names
	gen <- suppressWarnings(as.numeric(win_names))
	if (any(is.na(gen))) gen <- seq_len(W)
	data.frame(ASDSF = ASDSF, Generation = gen, row.names = NULL)
}

# Create short, deterministic IDs for clade keys (tip-index split strings)
shorten_clade_ids <- function(keys, prefix = "c_", n_hex = 16){
	if (!length(keys)) return(character(0L))
	hex <- vapply(keys, function(k) digest::digest(k, algo = "xxhash64"), "", USE.NAMES = FALSE)
	ids <- paste0(prefix, substr(hex, 1L, n_hex))
	# Ultra-rare: if truncation collides, append a sequence suffix to duplicates
	if (any(duplicated(ids))){
		dup_mask <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
		dup_vals <- ids[dup_mask]
		seq_in_grp <- ave(seq_along(dup_vals), dup_vals, FUN = function(z) seq_along(z))
		ids[which(dup_mask)] <- paste0(dup_vals, "-", seq_in_grp)
	}
	ids
}

# --- Edge-list native pipeline: clade frequencies and slide freq tables ---

# Compute clade frequencies for a window of edge-list trees.
# Returns a data.frame with columns: cladenames (character key), cladefreqs (numeric)
clade_freq_edges <- function(edge.list.window){
	stopifnot(is.list(edge.list.window))
	m <- length(edge.list.window)
	if (m == 0L){
		return(data.frame(cladenames = character(0), cladefreqs = numeric(0)))
	}
	# accumulate counts of canonical split keys across trees in the window
	acc <- new.env(parent = emptyenv())
	for (tr in edge.list.window){
		keys <- split_keys_from_edges(tr)
		if (!length(keys)) next
		for (k in keys){
			val <- if (exists(k, envir = acc, inherits = FALSE)) get(k, envir = acc) else 0L
			assign(k, val + 1L, envir = acc)
		}
	}
	nms <- ls(envir = acc, all.names = FALSE)
	if (!length(nms)){
		return(data.frame(cladenames = character(0), cladefreqs = numeric(0)))
	}
	cnts <- vapply(nms, function(k) get(k, envir = acc, inherits = FALSE), integer(1L))
	freq <- as.numeric(cnts) / m
	data.frame(cladenames = nms, cladefreqs = freq, stringsAsFactors = FALSE)
}

get_slide_freq_table_edges <- function(edge.chain, burnin = 0, window.size = 20, gens.per.tree = 1){
	stopifnot(is.list(edge.chain))
	if (burnin < 0L) burnin <- 0L
	# drop burnin once; keep only full windows
	edge.chain <- edge.chain[(burnin + 1L):length(edge.chain)]
	N <- length(edge.chain)
	if (!N) return(data.frame())
	n.windows <- as.integer(N / window.size)
	if (n.windows < 1L){
		# single (partial) window — keep behavior identical to previous version
		m <- N
		key_list <- vector("list", m)
		for (i in seq_len(m)) key_list[[i]] <- split_keys_from_edges(edge.chain[[i]])
		keys <- unlist(key_list, use.names = FALSE)
		tab <- table(keys)
		freq <- as.numeric(tab) / m
		vecs <- list(`1` = structure(freq, names = names(tab)))
		win_names <- prettyNum(((burnin + 1L) * gens.per.tree), sci = TRUE)
	} else {
		# precompute named frequency vectors per window using table(); no merges inside the loop
		vecs <- vector("list", n.windows)
		win_names <- prettyNum(seq((burnin + 1L), (burnin + n.windows * window.size), by = window.size) * gens.per.tree, sci = TRUE)
		# process each window without constructing intermediary data.frames
		for (w in seq_len(n.windows)){
			start <- (w - 1L) * window.size + 1L
			end <- start + window.size - 1L
			m <- window.size
			# collect per-tree unique keys then count with one table()
			key_list <- vector("list", m)
			kk <- 1L
			for (i in start:end){
				keys_i <- split_keys_from_edges(edge.chain[[i]])
				if (length(keys_i)){
					key_list[[kk]] <- keys_i
					kk <- kk + 1L
				}
			}
			if (kk == 1L){
				vecs[[w]] <- structure(numeric(0L), names = character(0L))
			} else {
				keys <- unlist(key_list[seq_len(kk - 1L)], use.names = FALSE)
				tab <- table(keys)
				vecs[[w]] <- structure(as.numeric(tab) / m, names = names(tab))
			}
		}
	}
	# union of clades across windows; build dense matrix (rows = clades, cols = windows)
	all_keys <- unique(unlist(lapply(vecs, names), use.names = FALSE))
	if (!length(all_keys)) return(data.frame())
	mat <- matrix(0, nrow = length(all_keys), ncol = length(vecs),
		dimnames = list(all_keys, win_names))
	for (w in seq_along(vecs)){
		v <- vecs[[w]]
		if (!length(v)) next
		mat[names(v), w] <- unname(v)
	}
	# attach short rownames + key map, compute sd/mean, order by sd desc
	short_ids <- shorten_clade_ids(rownames(mat))
	clade_map <- data.frame(id = short_ids, key = rownames(mat), stringsAsFactors = FALSE)
	df <- as.data.frame(mat, stringsAsFactors = FALSE, optional = TRUE)
	rownames(df) <- short_ids
	attr(df, "clade_key_map") <- clade_map
	df$sd <- apply(df, 1L, sd)
	df$mean <- rowMeans(as.matrix(df[, win_names, drop = FALSE]))
	df[order(df$sd, decreasing = TRUE), , drop = FALSE]
}


