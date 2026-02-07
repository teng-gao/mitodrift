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

#' Compute clade retention curve from tree
#' Returns data frame with confidence values and cumulative count of clades remaining
compute_retention_curve <- function(tree) {
    # Get confidence values from node labels
    conf_values <- as.numeric(tree$node.label)
    conf_values <- conf_values[!is.na(conf_values)]
    total_clades <- length(conf_values)

    # Sort descending to get unique thresholds
    conf_sorted <- sort(unique(conf_values), decreasing = TRUE)

    # For each unique confidence value, count how many clades remain at that threshold
    n_remaining <- vapply(conf_sorted, function(thresh) {
        sum(conf_values >= thresh)
    }, numeric(1))

    data.frame(
        confidence = conf_sorted,
        n_remaining = n_remaining,
        frac_remaining = n_remaining / total_clades
    )
}

#' Plot retention curves
#'
#' @param retention_data A single retention data frame, or a named list of them.
#'   When a list, names are used as legend labels.
#' @param sample_name Plot title.
#' @param cutoff Optional vertical line at a confidence threshold.
plot_retention_curve <- function(retention_data, sample_name = '', cutoff = NULL) {

    if (is.data.frame(retention_data)) {
        retention_data <- list(retention_data)
    }
    if (is.null(names(retention_data)) || any(names(retention_data) == "")) {
        base_name <- if (!is.null(sample_name) && nzchar(sample_name)) sample_name else "retention"
        if (length(retention_data) == 1L) {
            names(retention_data) <- base_name
        } else {
            names(retention_data) <- paste0(base_name, "_", seq_along(retention_data))
        }
    }
    df <- do.call(rbind, lapply(names(retention_data), function(nm) {
        d <- retention_data[[nm]]
        if (!all(c("confidence", "n_remaining") %in% colnames(d))) {
            stop("retention_data entries must contain columns: confidence, n_remaining")
        }
        d$group <- nm
        d
    }))

    single <- length(retention_data) == 1L

    if (single) {
        p <- ggplot(df, aes(x = confidence, y = n_remaining)) +
            geom_step(direction = 'hv', linewidth = 0.8, color = 'steelblue')
    } else {
        p <- ggplot(df, aes(x = confidence, y = n_remaining, color = group)) +
            geom_step(direction = 'hv', linewidth = 0.8)
    }

    if (!is.null(cutoff)) {
        cutoff <- as.numeric(cutoff)
        cutoff <- cutoff[!is.na(cutoff)]
        if (length(cutoff) == 0L) {
            stop("cutoff must contain at least one numeric value")
        }
        p <- p + geom_vline(xintercept = cutoff, color = 'red', linetype = 'dashed', linewidth = 0.6)
    }

    p <- p +
        scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
        labs(title = sample_name, x = 'Confidence threshold', y = 'Clades remaining', color = NULL) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            legend.position = if (single) 'none' else 'bottom',
            panel.grid.minor = element_blank()
        )

    return(p)
}

# Diagnostics utilities (precision/recall vs confidence)
#
# Dependencies (loaded by caller; do NOT library() here):
# - ape (for phylo handling; e.g., keep.tip)
# - ggplot2 (for plotting)
# - R/tree_evaluation.R (expects `variant_tree_jaccard_matrix()`)
# - mitodrift (expects `long_to_mat()`)
#
# These functions are extracted from `notebooks/Supplement/plot_all_prec_recall_vs_conf.R`
# and refactored to work directly from `phy_annot` and `mut_dat`.

compute_variant_pr_curve <- function(
    phy_annot,
    mut_dat,
    min_vaf = 0.01,
    min_cells = 2,
    j_thres = 0.66,
    n_points = 50
) {
    if (!inherits(phy_annot, "phylo")) stop("phy_annot must be a 'phylo' object")
    if (!exists("variant_tree_jaccard_matrix", mode = "function")) {
        stop("variant_tree_jaccard_matrix() not found. Source `R/tree_evaluation.R` before calling.")
    }
    if (!exists("long_to_mat", mode = "function")) {
        stop("long_to_mat() not found. `library(mitodrift)` (or otherwise load it) before calling.")
    }

    if (!all(c("variant", "cell", "a", "d") %in% colnames(mut_dat))) {
        stop("mut_dat must contain columns: variant, cell, a, d")
    }

    mut_dat <- mut_dat
    mut_dat$vaf <- mut_dat$a / mut_dat$d
    mut_dat$vaf[is.na(mut_dat$vaf) | is.infinite(mut_dat$vaf)] <- NA_real_

    # Provided by mitodrift: tools/mitodrift/R/utils.R
    vmat <- long_to_mat(mut_dat, "vaf")

    common_cells <- intersect(colnames(vmat), phy_annot$tip.label)
    common_cells <- sort(common_cells)
    if (length(common_cells) == 0L) stop("No overlapping cells between mut_dat and phy_annot tip labels.")

    vmat <- vmat[, common_cells, drop = FALSE]
    tr <- ape::keep.tip(phy_annot, common_cells)

    node_conf <- suppressWarnings(as.numeric(tr$node.label))
    node_conf[is.na(node_conf)] <- 0
    node_ids <- (length(tr$tip.label) + 1):(length(tr$tip.label) + tr$Nnode)
    names(node_conf) <- as.character(node_ids)

    J_full <- variant_tree_jaccard_matrix(
        vmat,
        tr,
        min_vaf = min_vaf,
        min_cells = min_cells
    )

    conf_all <- stats::setNames(node_conf, paste0("Node", names(node_conf)))

    df <- compute_pr_from_jaccard_mat(J_full, conf_all, j_thres = j_thres, n_points = n_points)

    var_max_j <- apply(J_full, 1, max)
    df$fraction_well_matched <- mean(var_max_j >= j_thres)
    df
}

#' Compute precision/recall/F1 curve from a truth tree and an estimated tree
#'
#' @param phy_tru A rooted \code{phylo} object representing the ground-truth tree.
#' @param phy_est A rooted \code{phylo} object with clade confidence in \code{node.label}.
#' @param j_thres Jaccard threshold for declaring a match (default 0.5).
#' @param n_points Number of confidence cutoff points to evaluate (default 50).
#'
#' @return Data frame with columns: conf_cutoff, n_pred, n_tp, n_tru,
#'   n_recover, precision, recall, f1.
compute_tree_pr_curve <- function(phy_tru, phy_est, j_thres = 0.5, n_points = 50) {
    if (!inherits(phy_tru, "phylo")) stop("phy_tru must be a 'phylo' object")
    if (!inherits(phy_est, "phylo")) stop("phy_est must be a 'phylo' object")

    common_tips <- sort(intersect(phy_tru$tip.label, phy_est$tip.label))
    if (length(common_tips) == 0L) stop("Trees share no tip labels.")

    phy_tru <- ape::keep.tip(phy_tru, common_tips)
    phy_est <- ape::keep.tip(phy_est, common_tips)

    J <- clade_jaccard_matrix(phy_tru, phy_est, add_node_names = FALSE)

    # Map column node IDs to confidence values
    n_tips <- length(phy_est$tip.label)
    col_ids <- as.integer(colnames(J))
    conf_all <- stats::setNames(
        suppressWarnings(as.numeric(phy_est$node.label[col_ids - n_tips])),
        colnames(J)
    )
    conf_all[is.na(conf_all)] <- 0

    compute_pr_from_jaccard_mat(J, conf_all, j_thres = j_thres, n_points = n_points)
}

#' Compute precision/recall/F1 curve from a Jaccard matrix and clade confidences
#'
#' @param J Matrix of Jaccard scores (rows = truth clades, cols = predicted clades).
#' @param conf Named numeric vector of clade confidence values keyed by node ID.
#'   Names must match column names of \code{J}.
#' @param j_thres Jaccard threshold for declaring a match (default 0.66).
#' @param n_points Number of confidence cutoff points to evaluate (default 50).
#'
#' @return Data frame with columns: conf_cutoff, n_pred, n_tp, n_tru,
#'   n_recover, precision, recall, f1.
compute_pr_from_jaccard_mat <- function(J, conf, j_thres = 0.66, n_points = 50) {

    conf_cols <- conf[colnames(J)]
    conf_cols[is.na(conf_cols)] <- 0

    # Per predicted clade: TP if best-matching truth clade >= j_thres
    col_supported <- apply(J, 2, max) >= j_thres

    # Per truth clade: max confidence among predicted clades that match it
    row_support_conf <- apply(J, 1, function(x) {
        idx <- which(x >= j_thres)
        if (length(idx) == 0L) return(-Inf)
        max(conf_cols[idx], na.rm = TRUE)
    })

    # Sweep confidence cutoffs
    probs <- seq(0, 1, length.out = n_points)
    conf_cutoffs <- as.numeric(stats::quantile(conf_cols, probs = probs, na.rm = TRUE))
    conf_cutoffs[1] <- 0
    conf_cutoffs[length(conf_cutoffs)] <- 1
    conf_cutoffs <- pmin(pmax(conf_cutoffs, 0), 1)

    n_tru <- nrow(J)
    score_mat <- sapply(conf_cutoffs, function(ct) {
        keep <- conf_cols >= ct
        n_pred <- sum(keep)
        n_tp <- sum(col_supported & keep)
        n_recover <- sum(row_support_conf >= ct)
        precision <- if (n_pred == 0L) NA_real_ else n_tp / n_pred
        recall <- if (n_tru == 0L) NA_real_ else n_recover / n_tru
        c(n_pred = n_pred, n_tp = n_tp, n_tru = n_tru,
          n_recover = n_recover, precision = precision, recall = recall)
    })

    df <- data.frame(
        conf_cutoff = conf_cutoffs,
        n_pred = score_mat["n_pred", ],
        n_tp = score_mat["n_tp", ],
        n_tru = score_mat["n_tru", ],
        n_recover = score_mat["n_recover", ],
        precision = score_mat["precision", ],
        recall = score_mat["recall", ]
    )
    denom <- df$precision + df$recall
    df$f1 <- ifelse(is.na(denom) | denom == 0, NA_real_,
                    2 * df$precision * df$recall / denom)
    df
}

plot_prec_recall_vs_conf <- function(
    pr_df,
    sample_name = NULL,
    cutoff = NULL,
    x_axis_scale = c("log10", "linear"),
    log10_min = 1e-6,
    j_thres = NULL
) {
    x_axis_scale <- match.arg(x_axis_scale)
    if (!all(c("conf_cutoff", "precision", "recall", "f1") %in% colnames(pr_df))) {
        stop("pr_df must contain columns: conf_cutoff, precision, recall, f1")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is required for plot_prec_recall_vs_conf(); please install/load it.")
    }

    fraction <- pr_df$fraction_well_matched
    fraction <- if (length(fraction) == 0L) NA_real_ else fraction[1]
    
    if (!is.null(j_thres)) {
        ann <- if (is.na(fraction)) NA_character_ else sprintf("%.1f%% vars maxJ >= %.2f", 100 * fraction, j_thres)
    } else {
        ann = ''
    }

    df_long <- rbind(
        data.frame(conf_cutoff = pr_df$conf_cutoff, metric = "Precision", value = pr_df$precision),
        data.frame(conf_cutoff = pr_df$conf_cutoff, metric = "Recall", value = pr_df$recall),
        data.frame(conf_cutoff = pr_df$conf_cutoff, metric = "F1", value = pr_df$f1)
    )

    x_lab <- "Confidence cutoff"
    vline_x <- cutoff
    x_var <- "conf_cutoff"
    if (x_axis_scale == "log10") {
        df_long$conf_cutoff_plot <- pmax(df_long$conf_cutoff, log10_min)
        x_var <- "conf_cutoff_plot"
        x_lab <- "Confidence cutoff (log10)"
        if (!is.null(vline_x)) vline_x <- max(vline_x, log10_min)
    }

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[x_var]], y = value, color = metric)) +
        ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
        ggplot2::scale_color_manual(values = c("Precision" = "#1F77B4", "Recall" = "#D62728", "F1" = "#2CA02C")) +
        ggplot2::labs(
            title = sample_name,
            x = x_lab,
            y = NULL
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
            axis.title = ggplot2::element_text(size = 8),
            axis.text = ggplot2::element_text(size = 7),
            legend.position = "bottom",
            panel.grid.minor = ggplot2::element_blank()
        )

    if (!is.null(vline_x)) {
        p <- p + ggplot2::geom_vline(
            xintercept = vline_x,
            color = "grey40",
            linetype = "dashed",
            linewidth = 0.6
        )
    }

    if (!is.na(ann)) {
        x_anno <- min(df_long[[x_var]], na.rm = TRUE)
        y_anno <- min(df_long$value, na.rm = TRUE)
        p <- p + ggplot2::annotate(
            "text",
            x = x_anno,
            y = y_anno,
            label = ann,
            hjust = 0,
            vjust = 0,
            size = 2.6
        )
    }

    if (x_axis_scale == "log10") {
        p <- p + ggplot2::scale_x_continuous(trans = "log10")
    }

    p
}
