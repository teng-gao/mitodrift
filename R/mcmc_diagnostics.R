## -----------------------------------------------------------------------------
#' Compute ASDSF across chains using clades from a target tree
#'
#' @param phy_target A rooted `phylo` object defining the reference clades. Assumed to be in postorder.
#' @param edge_list_chains List of chains, each a list of edge matrices (2-column integer matrices).
#' @param rooted Logical; treat trees as rooted when matching clades. Default `TRUE`.
#' @param normalize Logical; pass through to `prop_clades_par` (default `TRUE`).
#' @param min_freq Minimum clade frequency threshold used when averaging SDs. Default `0.1`.
#'
#' @return A list with elements `asdsf`, `per_clade_sd`, `keep_mask`, and `freq_matrix`.
#' @keywords internal
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


