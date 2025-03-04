#' @import dplyr
#' @import tidygraph
#' @import stringr
#' @importFrom igraph vcount ecount E V V<- E<- 
#' @importFrom phangorn upgma 
#' @importFrom ape root drop.tip nj
#' @importFrom parallelDist parDist
#' @importFrom stats na.omit reorder setNames
#' @useDynLib scistreer
NULL

#' Run the scistree workflow
#' @param P matrix Genotype probability matrix (cell x mutation). Each entry is a probability (0-1) that the cell harbors the mutation
#' @param init character Initialization strategy; UPGMA or NJ
#' @param max_iter integer Maximum number of iterations
#' @param eps numeric Tolerance threshold in likelihood difference for stopping
#' @param verbose logical Verbosity
#' @param ncores integer Number of cores to use
#' @return phylo A maximum-likelihood phylogeny
#' @examples
#' tree_small = run_scistree(P_small)
#' @export
run_scistree = function(P, init = 'UPGMA', ncores = 1, max_iter = 100, eps = 0.01, verbose = TRUE) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)

    if (!is.matrix(P)) {
        stop('P needs to be a matrix')
    }

    # contruct initial tree
    dist_mat = parDist(rbind(P, 'outgroup' = 0), threads = ncores)

    if (init == 'UPGMA') {

        tree_init = upgma(dist_mat) %>%
            root(outgroup = 'outgroup') %>%
            drop.tip('outgroup') %>%
            reorder(order = 'postorder')
            
    } else if (init == 'NJ') {

        tree_init = nj(dist_mat) %>%
            root(outgroup = 'outgroup') %>%
            drop.tip('outgroup') %>%
            reorder(order = 'postorder')

    } else {
        stop("init has be one of 'UPGMA', 'NJ'")
    }

    tree_list = perform_nni(tree_init, P, ncores = ncores, max_iter = max_iter, eps = eps, verbose = verbose)

    tree_final = tree_list[[length(tree_list)]]

    # get the branch lengths
    tree_final$edge.length = rep(1, nrow(tree_final$edge))

    return(tree_final)
    
}


#' Maximum likelihood tree search via NNI
#' @param tree_init phylo Intial tree
#' @param P matrix Genotype probability matrix
#' @param max_iter integer Maximum number of iterations
#' @param eps numeric Tolerance threshold in likelihood difference for stopping
#' @param verbose logical Verbosity
#' @param ncores integer Number of cores to use
#' @return multiPhylo List of trees corresponding to the rearrangement steps
#' @examples
#' tree_list = perform_nni(tree_upgma, P_small)
#' @export 
perform_nni = function(tree_init, P, max_iter = 100, eps = 0.01, ncores = 1, verbose = TRUE) {

    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
        
    converge = FALSE
    
    i = 1
    max_current = -Inf
    tree_current = tree_init
    tree_list = list()
    tree_list[[1]] = tree_current
    
    while (!converge & i <= max_iter) {
        
        i = i + 1
        
        ptm = proc.time()
        
        RcppParallel::setThreadOptions(numThreads = ncores)
        
        scores = nni_cpp_parallel(tree_current, P)
        
        if (max(scores) > max_current + eps) {
            max_id = which.max(scores)
            if (max_id %% 2 == 0) {pair_id = 2} else {pair_id = 1}
            tree_current$edge = matrix(nnin_cpp(tree_current$edge, ceiling(max_id/2))[[pair_id]], ncol = 2)
            tree_list[[i]] = tree_current
            tree_list[[i]]$likelihood = max_current = max(scores)
            converge = FALSE
        } else {
            converge = TRUE
        }
        
        runtime = proc.time() - ptm
        
        if (verbose) {
            msg = paste0('Iter ', i, ' ', max_current, ', ', signif(unname(runtime[3]),2), 's')
            message(msg)
        }
    }
    
    class(tree_list) = 'multiPhylo'
    
    return(tree_list)
}

#' Score a tree based on maximum likelihood
#' @param tree phylo object
#' @param P genotype probability matrix
#' @param get_l_matrix whether to compute the whole likelihood matrix
#' @return list Likelihood scores of a tree
#' @examples
#' tree_likelihood = score_tree(tree_upgma, P_small)$l_tree
#' @export
score_tree = function(tree, P, get_l_matrix = FALSE) {
 
    tree = reorder(tree, order = 'postorder')
        
    n = nrow(P)
    m = ncol(P)

    logQ = matrix(nrow = tree$Nnode * 2 + 1, ncol = m)

    logP_0 = log(1-P)
    logP_1 = log(P)
    
    node_order = c(tree$edge[,2], n+1)
    node_order = node_order[node_order > n]
    
    logQ[1:n,] = logP_1 - logP_0

    children_dict = allChildrenCPP(tree$edge)

    logQ = CgetQ(logQ, children_dict, node_order)

    if (get_l_matrix) {
        l_matrix = sweep(logQ, 2, colSums(logP_0), FUN = '+')
        l_tree = sum(apply(l_matrix, 2, max))
    } else {
        l_matrix = NULL
        l_tree = sum(apply(logQ, 2, max)) + sum(logP_0)
    }
    
    return(list('l_tree' = l_tree, 'logQ' = logQ, 'l_matrix' = l_matrix))
    
}


#' Find maximum lilkelihood assignment of mutations on a tree
#' @param tree phylo Single-cell phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @return tbl_graph A single-cell phylogeny with mutation placements
#' @examples
#' gtree_small = annotate_tree(tree_small, P_small)
#' @export
annotate_tree = function(tree, P) {
   
    sites = colnames(P)
    n = nrow(P)
    tree_stats = score_tree(tree, P, get_l_matrix = TRUE)

    l_matrix = as.data.frame(tree_stats$l_matrix)

    colnames(l_matrix) = sites
    rownames(l_matrix) = c(tree$tip.label, paste0('Node', 1:tree$Nnode))

    # mutation assignment on nodes
    mut_nodes = data.frame(
            site = sites,
            node_phylo = apply(l_matrix, 2, which.max),
            l = apply(l_matrix, 2, max)
        ) %>%
        mutate(
            name = ifelse(node_phylo <= n, tree$tip.label[node_phylo], paste0('Node', node_phylo - n))
        ) %>%
        group_by(name) %>%
        summarise(
            site = paste0(sort(site), collapse = ','),
            n_mut = n(),
            l = sum(l),
            .groups = 'drop'
        )

    gtree = tree %>%
        ladderize() %>%
        as_tbl_graph() %>%
        mutate(
            leaf = node_is_leaf(),
            root = node_is_root(),
            depth = bfs_dist(root = 1),
            id = row_number()
        )

    # leaf annotation for edges
    gtree = gtree %>%
        activate(edges) %>%
        select(-any_of(c('leaf'))) %>%
        left_join(
            gtree %>%
                activate(nodes) %>%
                data.frame() %>%
                select(id, leaf),
            by = c('to' = 'id')
        )

    # annotate the tree
    gtree = mut_to_tree(gtree, mut_nodes)

    return(gtree)
}

#' Transfer mutation assignment onto a single-cell phylogeny
#'
#' @param gtree tbl_graph The single-cell phylogeny
#' @param mut_nodes dataframe Mutation placements
#' @return tbl_graph A single-cell phylogeny with mutation placements
#' @examples
#' gtree_small = mut_to_tree(gtree_small, mut_nodes_small)
#' @export
mut_to_tree = function(gtree, mut_nodes) {
   
    # transfer mutation to tree
    gtree = gtree %>%
        activate(nodes) %>%
        select(-any_of(c('n_mut', 'l', 'site', 'clone'))) %>%
        left_join(
            mut_nodes %>%
                mutate(n_mut = unlist(lapply(str_split(site, ','), length))) %>%
                select(name, n_mut, site),
            by = 'name'
        ) %>%
        mutate(n_mut = ifelse(is.na(n_mut), 0, n_mut))

    # get branch length
    gtree = gtree %>% 
        activate(edges) %>%
        select(-any_of(c('length'))) %>%
        left_join(
            gtree %>%
                activate(nodes) %>%
                data.frame() %>%
                select(id, length = n_mut),
            by = c('to' = 'id')
        ) %>%
        mutate(length = ifelse(leaf, pmax(length, 0.2), length))

    # label genotype on nodes
    node_to_mut = gtree %>% activate(nodes) %>% data.frame() %>% {setNames(.$site, .$id)}

    gtree = gtree %>%
        activate(nodes) %>%
        mutate(GT = unlist(
            map_bfs(node_is_root(),
            .f = function(path, ...) { paste0(na.omit(node_to_mut[path$node]), collapse = ',') })
            ),
            last_mut = unlist(
                map_bfs(node_is_root(),
                .f = function(path, ...) { 
                    past_muts = na.omit(node_to_mut[path$node])
                    if (length(past_muts) > 0) {
                        return(past_muts[length(past_muts)])
                    } else {
                        return('')
                    }
                })
            )
        ) %>%
        mutate(GT = ifelse(GT == '' & !is.na(site), site, GT))

    # preserve the clone ids
    if ('GT' %in% colnames(mut_nodes)) {
        gtree = gtree %>% activate(nodes) %>%
            left_join(
                mut_nodes %>% select(GT, clone),
                by = 'GT'
            )
    }
    
    return(gtree)
}

#' Convert a single-cell phylogeny with mutation placements into a mutation graph
#'
#' @param gtree tbl_graph The single-cell phylogeny
#' @return igraph Mutation graph
#' @examples
#' mut_graph = get_mut_graph(gtree_small)
#' @export
get_mut_graph = function(gtree) {

    mut_nodes = gtree %>% 
        activate(nodes) %>%
        as.data.frame() %>%
        filter(!is.na(site)) %>%
        distinct(name, site)

    G = gtree %>%
        activate(nodes) %>%
        arrange(last_mut) %>%
        convert(to_contracted, last_mut) %>%
        mutate(label = last_mut, id = 1:n()) %>%
        as.igraph

    G = label_edges(G)

    V(G)$node = G %>%
        igraph::as_data_frame('vertices') %>%
        left_join(
            mut_nodes %>% dplyr::rename(node = name),
            by = c('label' = 'site')
        ) %>%
        pull(node)

    G = G %>% transfer_links()

    return(G)

}


#' Annotate the direct upstream or downstream mutations on the edges
#'
#' @param G igraph Mutation graph
#' @return igraph Mutation graph 
#' @keywords internal
label_edges = function(G) {
    
    edge_df = G %>% igraph::as_data_frame('edges') %>%
        left_join(
            G %>% igraph::as_data_frame('vertices') %>% select(from_label = label, id),
            by = c('from' = 'id')
        ) %>%
        left_join(
            G %>% igraph::as_data_frame('vertices') %>% select(to_label = label, id),
            by = c('to' = 'id')
        ) %>%
        mutate(label = paste0(from_label, '->', to_label))
    
    E(G)$label = edge_df$label
    E(G)$from_label = edge_df$from_label
    E(G)$to_label = edge_df$to_label
    
    return(G)
}


#' Annotate the direct upstream or downstream node on the edges
#'
#' @param G igraph Mutation graph
#' @return igraph Mutation graph 
#' @keywords internal
transfer_links = function(G) {
    
    edge_df = G %>% igraph::as_data_frame('edges') %>%
            left_join(
                G %>% igraph::as_data_frame('vertices') %>% select(from_node = node, id),
                by = c('from' = 'id')
            ) %>%
            left_join(
                G %>% igraph::as_data_frame('vertices') %>% select(to_node = node, id),
                by = c('to' = 'id')
            )

    E(G)$from_node = edge_df$from_node
    E(G)$to_node = edge_df$to_node
    
    return(G)
}

#' Convert the phylogeny from tidygraph to phylo object
#' modified from R package alakazam, converts a tbl_graph to a phylo object
#'
#' @param graph tbl_graph The single-cell phylogeny
#' @return phylo The single-cell phylogeny
#' @examples
#' tree_small = to_phylo(annotate_tree(tree_small, P_small))
#' @export
to_phylo = function(graph) {
    
    df  <- igraph::as_data_frame(graph)
    node_counts <- table(c(df$to,df$from))
    tips <- names(node_counts)[node_counts == 1]
    nodes <- names(node_counts)[node_counts > 1]
    attr <- igraph::vertex_attr(graph)

    tipn <- 1:length(tips)
    names(tipn) <- tips
    noden <- (length(tips)+1):(length(tips)+length(nodes))
    names(noden) <- nodes
    renumber <- c(tipn,noden)

    df$from <- as.numeric(renumber[df$from])
    df$to <- as.numeric(renumber[df$to])    

    phylo <- list()
    phylo$edge <- matrix(cbind(df$from,df$to),ncol=2)
    phylo$edge.length <- as.numeric(df$length)
    phylo$tip.label <- tips
    phylo$Nnode <- length(nodes)
    class(phylo) <- "phylo"

    nnodes <- length(renumber)
    phylo$nodes <- lapply(1:nnodes,function(x){
        n <- list()
        n$id <- names(renumber[renumber == x])
        n
    })
    
    n_mut_root = graph %>% activate(nodes) %>% filter(node_is_root()) %>% pull(n_mut)
    
    phylo$root.edge = n_mut_root
    
    return(phylo)
}


#' From ape; will remove once new ape version is released
#' https://github.com/emmanuelparadis/ape/issues/54
#' @param phy phylo The phylogeny
#' @param right logical Whether ladderize to the right
#' @examples
#' tree_small = ladderize(tree_small)
#' @export
ladderize <- function(phy, right = TRUE) {
    desc_fun <- function(x) {
        parent <- x[, 1]
        children <- x[, 2]
        res <- vector("list", max(x))
        for (i in seq_along(parent)) res[[parent[i]]] <- c(res[[parent[i]]], children[i])
        return(res)
    }

    if(!is.null(phy$edge.length)){
        el <- numeric(max(phy$edge))
        el[phy$edge[, 2]] <- phy$edge.length
    }

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    nb.edge <- dim(phy$edge)[1]

    phy <- reorder(phy, "postorder")
    N <- node_depth(as.integer(nb.tip), as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
            as.integer(nb.edge), double(nb.tip + nb.node), 1L)

    ii <- order(x <- phy$edge[,1], y <- N[phy$edge[,2]], decreasing = right)
    desc <- desc_fun(phy$edge[ii,])

    tmp <- integer(nb.node)
    new_anc <- integer(nb.node)
    new_anc[1] <- tmp[1] <- nb.tip + 1L
    k <- nb.node
    pos <- 1L

    while(pos > 0L && k > 0){
        current <- tmp[pos]
        new_anc[k] <- current
        k <- k - 1L
        dc <- desc[[current]]
        ind <- (dc > nb.tip)
        if(any(ind)){
            l <- sum(ind)
            tmp[pos -1L + seq_len(l)] <-  dc[ind]
            pos <- pos + l - 1L
        }
        else pos <- pos - 1L
    }
    edge <- cbind(rep(new_anc, lengths(desc[new_anc])), unlist(desc[new_anc]))
    phy$edge <- edge
    if(!is.null(phy$edge.length)) phy$edge.length <- el[edge[,2]]
    attr(phy, "order") <- "postorder"
    phy <- reorder(phy, "cladewise")
    phy
}
