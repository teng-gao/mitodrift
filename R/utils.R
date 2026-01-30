
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

#' Collapse weak clades below a confidence threshold.
#'
#' Given a fully binary phylogeny with per-node confidence scores stored in
#' `tree$node.label`, this helper collapses every internal node whose confidence
#' falls below `conf`, effectively pruning low-support clades while retaining
#' higher-confidence structure.
#'
#' @param tree A binary `phylo` object whose `node.label` vector stores
#'   posterior/confidence values for internal nodes.
#' @param conf Numeric threshold in $[0,1]$. Internal nodes with confidence below
#'   this value are collapsed.
#' @return A renumbered `phylo` object with low-confidence nodes collapsed.
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

#' Collapse branches with high expected mis-assignments.
#'
#' Computes the expected number of incorrectly assigned tips for each internal
#' branch (based on $(1 - p) \times$ clade size, where $p$ is the branch
#' confidence) and collapses branches whose expectation exceeds a tolerance
#' derived from the total number of tips.
#'
#' @param tree A binary `phylo` object with internal node confidence scores in
#'   `tree$node.label`.
#' @param tol Numeric tolerance expressed as a fraction of total tips; the
#'   threshold is `length(tree$tip.label) * tol` expected errors.
#' @return A renumbered `phylo` object with overly uncertain branches collapsed.
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

#' Assign clone IDs to a tree allowing small polytomies.
#'
#' @param tree Rooted `phylo` object with bifurcating or polytomous structure.
#' @param k Positive scalar; maximum clade size that collapses into one clone.
#' @param paraphyletic Logical; if `TRUE`, leftover unassigned tips below a node
#'   are grouped into one clone, otherwise each becomes a singleton.
#' @param return_df Logical; if `TRUE`, returns a data frame annotated per tip,
#'   otherwise an integer vector of clone assignments.
#'
#' @return Either a data frame with columns `cell`, `clade`, `annot`, `size`,
#'   `frac` or an integer vector of clone IDs indexed by tip order.
#'
#' @details Finds the root (node never used as a child), recurses over internal
#' nodes, and assigns clone IDs whenever a clade has `<= k` tips. Tips directly
#' attached to the root are always singleton clones. Ensures all tips receive
#' an assignment, warning if late singletons are required.
#' @export
assign_clones_polytomy <- function(tree, k = Inf, paraphyletic = FALSE, return_df = TRUE) {

	if (!inherits(tree, "phylo")) {
		stop("tree must be a 'phylo' object")
	}
	if (!is.numeric(k) || length(k) != 1L || k <= 0) {
		stop("k must be a positive scalar")
	}

	Ntip <- length(tree$tip.label)
	if (Ntip == 0L) {
		return(integer(0L))
	}

	# Find root as the node that never appears as a child
	all_nodes <- unique(as.vector(tree$edge))
	child_nodes <- unique(tree$edge[, 2])
	root_node <- setdiff(all_nodes, child_nodes)
	if (length(root_node) != 1L) {
		stop("Could not uniquely identify root node")
	}
	root_node <- root_node[1L]

	total_nodes <- Ntip + tree$Nnode

	# children_list: parent -> vector of children node indices
	children_list <- split(tree$edge[, 2], tree$edge[, 1])

	# Precompute descendant tips for each node (memoized)
	subtree_tips <- vector("list", total_nodes)
	get_subtree_tips <- function(nd) {
		if (nd <= Ntip) return(nd)
		if (!is.null(subtree_tips[[nd]])) return(subtree_tips[[nd]])

		ch <- children_list[[as.character(nd)]]
		if (is.null(ch) || length(ch) == 0L) {
			subtree_tips[[nd]] <<- integer(0L)
			return(integer(0L))
		}
		tips_nd <- integer(0L)
		for (cnd in ch) {
			if (cnd <= Ntip) {
				tips_nd <- c(tips_nd, cnd)
			} else {
				tips_nd <- c(tips_nd, get_subtree_tips(cnd))
			}
		}
		tips_nd <- unique(tips_nd)
		subtree_tips[[nd]] <<- tips_nd
		tips_nd
	}

	# clone_assign[i] = clone ID for tip i (0 = unassigned)
	clone_assign <- integer(Ntip)
	clone_node_assign <- integer(Ntip)  # track which node defines each clone
	next_clone <- 1L

	assign_tips <- function(tips, node_id) {
		tips <- as.integer(tips)
		if (length(tips) == 0L) return(invisible(NULL))
		# only assign still-unassigned tips
		tips <- tips[clone_assign[tips] == 0L]
		if (length(tips) == 0L) return(invisible(NULL))
		clone_assign[tips] <<- next_clone
		clone_node_assign[tips] <<- node_id
		next_clone <<- next_clone + 1L
		invisible(NULL)
	}

	# Recursive processing for internal nodes below the root
	process_node <- function(nd) {
		# nd is guaranteed > Ntip here
		node_tips <- get_subtree_tips(nd)

		# If this whole clade is small enough, make it a single clone
		if (length(node_tips) <= k) {
			assign_tips(node_tips, nd)
			return(invisible(NULL))
		}

		# Otherwise, recurse into internal children first
		ch <- children_list[[as.character(nd)]]
		if (!is.null(ch) && length(ch) > 0L) {
			for (cnd in ch) {
				if (cnd > Ntip) {
					process_node(cnd)
				}
			}
		}

		# After children have formed their clones,
		# group any leftover tips under this node
		remaining <- node_tips[clone_assign[node_tips] == 0L]
		if (length(remaining) > 0L) {
			if (paraphyletic) {
				# group all remaining paraphyletic tips into a single clone
				assign_tips(remaining, nd)
			} else {
				# assign each remaining tip as its own singleton clone
				for (tip in remaining) {
					assign_tips(tip, tip)
				}
			}
		}

		invisible(NULL)
	}

	# --- 1) Handle root-level children ---
	root_children <- children_list[[as.character(root_node)]]

	# a) Tips directly from root are always their own singleton clones
	for (child in root_children) {
		if (child <= Ntip) {
			assign_tips(child, child)
		}
	}

	# b) Internal children are processed with the recursive logic
	for (child in root_children) {
		if (child > Ntip) {
			process_node(child)
		}
	}

	# Sanity: everything should be assigned
	if (any(clone_assign == 0L)) {
		warning("Some tips remain unassigned; they will be given unique singleton clones.")
		unassigned <- which(clone_assign == 0L)
		for (i in unassigned) {
			assign_tips(i, i)
		}
	}

	names(clone_assign) <- tree$tip.label
	names(clone_node_assign) <- tree$tip.label

	if (return_df) {
		clade_annot = data.frame(clade = clone_assign, clade_node = clone_node_assign) %>%
			tibble::rownames_to_column("cell") %>%
			mutate(annot = as.character(clade)) %>%
			group_by(clade) %>%
			mutate(size = n()) %>%
			ungroup() %>%
			mutate(frac = size / n()) %>%
			arrange(-frac) %>%
            mutate(clade = as.integer(factor(clade, unique(clade)))) %>%
            mutate(clade = as.character(clade)) %>%
            mutate(annot = clade)
		return(clade_annot)
	} else {
		return(list(
			clade = clone_assign,
			clade_node = clone_node_assign
		))
	}
}