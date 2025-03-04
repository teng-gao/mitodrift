#' @import ggplot2
#' @importFrom ggtree ggtree geom_rootedge
#' @import patchwork
NULL

#' Plot phylogeny and mutation heatmap
#' @param tree phylo or tbl_graph Phylogeny
#' @param P matrix Genotype probability matrix
#' @param branch_width numeric Branch width
#' @param root_edge logical Whether to plot root edge
#' @return ggplot Plot visualizing the single-cell phylogeny and mutation probability heatmap
#' @examples
#' p = plot_phylo_heatmap(tree_small, P_small)
#' @export
plot_phylo_heatmap = function(tree, P, branch_width = 0.5, root_edge = TRUE) {

    if ('tbl_graph' %in% class(tree)) {
        tree = to_phylo(tree)
    }

    # plot phylogeny 
    p_tree = tree %>%
        ggtree(ladderize = TRUE, size = branch_width) +
        theme(
            plot.margin = margin(0,1,0,0, unit = 'mm'),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            # axis.text.y = element_text(size = 5)
            axis.text.y = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent", color = NA)
        ) +
        guides(color = 'none')

    if (root_edge) {
        p_tree = p_tree + geom_rootedge(size = branch_width)
    }

    # order the cells
    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

    p_heatmap = as.data.frame(P) %>%
        mutate(cell = rownames(P)) %>%
        mutate(cell = factor(cell, cell_order)) %>%
        reshape2::melt(id.var = 'cell', variable.name = 'mut', value.name = 'p') %>%
        ggplot(
            aes(x = mut, y = cell, fill = p)
        ) +
        geom_tile() +
        theme_void() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.ticks.x = element_line()
        )

    (p_tree | p_heatmap) + plot_layout(widths = c(1,2))

}
