#' @export 
plot_phylo_heatmap2 = function(gtree, df_var, branch_width = 0.25, root_edge = TRUE, dot_size = 1, ylim = NULL,
    clade_annot = NULL,
    title = NULL, auc = FALSE, clone_bar = FALSE, label_site = FALSE, cell_annot = NULL, tip_lab = FALSE, node_lab = FALSE,
    het_max = 0.1, conf_min = 0.5, conf_label = FALSE, branch_length = TRUE, node_conf = FALSE, annot_scale = NULL, annot_legend = FALSE, label_group = FALSE,
    annot_legend_title = '', text_size = 3, label_size = 1, mut = NULL, post_max = FALSE, mark_low_cov = FALSE, facet_by_group = FALSE, flip = FALSE) {

    if (inherits(gtree, 'tbl_graph')) {
        phylo = to_phylo_reorder(gtree)
    } else {
        phylo = gtree
        node_conf = FALSE
    }
    
    if (!branch_length) {
        phylo$edge.length = NULL
    }

    p_tree = phylo %>%
        ggtree(ladderize = TRUE, linewidth = branch_width) + 
        theme_bw() +
        theme(
            plot.margin = margin(0, 0.1, 0, 0, unit = "mm"), 
            axis.title.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x = element_blank(),
            axis.text.y = element_text(), 
            axis.line.y = element_line(),
            axis.ticks.y = element_line(), 
            panel.background = element_rect(fill = "transparent", colour = NA), 
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid = element_blank()
        ) + 
        coord_flip() +
        scale_x_reverse() +
        scale_y_continuous(expand = expansion(add = 1)) +
        ggtitle(title)

    if (flip) {
        p_tree = p_tree + scale_y_reverse(expand = expansion(add = 1))
    }

    # plot mutation VAF
    if (!is.null(mut)) {

        dat = gtree %>% activate(nodes) %>%
            as.data.frame() %>% select(all_of(c('name', p_v = mut)))

        p_tree = p_tree %<+% 
            dat + 
            ggraph::geom_node_point(aes(color = p_v), size = dot_size) +
            scale_color_gradient(limits = c(0,het_max), oob = scales::oob_squish)

    }

    if (node_conf) {

        dat = gtree %>% activate(nodes) %>%
            mutate(isRoot = node_is_root()) %>%
            as.data.frame() %>% 
            select(any_of(c('name', 'isRoot', 'conf')))

        if ('conf' %in% colnames(dat)) {
            p_tree = p_tree %<+% 
                dat + 
                geom_nodepoint(aes(fill = conf, subset = !isTip & !isRoot, x = branch), size = dot_size, pch = 22, stroke = 0) +
                scale_fill_gradient(low = 'white', high = 'firebrick', limits = c(conf_min,1), oob = scales::oob_squish)

            if (conf_label) {
                p_tree = p_tree + 
                    geom_text2(
                        aes(label = round(conf, 2), x = branch, subset = !isTip & !isRoot),
                        size = text_size, hjust = 0, vjust = 0.25
                    )
            }
        }

    }

    if (tip_lab) {
        p_tree = p_tree + geom_tiplab(size = label_size, vjust = 0.5, hjust = 1.2, angle = 90) +
            scale_x_reverse(expand = expansion(mult = c(0.15, 0.05)))
    }

    if (node_lab) {
        p_tree = p_tree + 
            geom_text2(
                aes(label = label, subset = !isTip, x = branch),
                size = label_size, vjust = 0.5, hjust = 0.5
            )
    }

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% 
        pull(label)
    
    df_var = df_var %>% filter(cell %in% cell_order) %>%
        mutate(cell = factor(as.integer(factor(cell, cell_order)), 1:length(cell_order)))

    if (!'vaf' %in% colnames(df_var)) {
        df_var = df_var %>% mutate(vaf = a/d)
    }

    p_heatmap = df_var %>% 
        filter(vaf > 0) %>%
        ggplot(
            aes(x = cell, y = variant, fill = vaf)
        ) +
        geom_tile(width = 0.8) + 
        theme_bw() + 
        theme( 
            plot.margin = margin(0, 1, 0, 0, unit = "mm"), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_text(size = text_size),
            panel.background = element_rect(fill = "transparent", colour = NA), 
            plot.background = element_rect(fill = "transparent"),
            # panel.grid.minor.y = element_blank(),
            # panel.grid.major.y = element_blank(),
            panel.grid = element_blank()
        ) +
        scale_x_discrete(expand = expansion(add = 1), drop = F) +
        scale_y_discrete(expand = expansion(mult = 0.01)) +
        scale_fill_gradient(low = 'white', high = 'red', limits = c(0,het_max), oob = scales::oob_squish) +
        guides(fill = guide_colorbar(title = 'VAF'))

    if (facet_by_group) {
        p_heatmap = p_heatmap + facet_grid(group~., scales = 'free_y', space = 'free_y') +
            theme(panel.spacing.y = unit(0,'mm'))
    }
    
    if (mark_low_cov) {
        p_heatmap = p_heatmap + geom_point(
            data = df_var %>% filter(d < 10),
            pch = 4,
            size = dot_size
        )
    }

    if (!is.null(clade_annot)) {
        p_clade = clade_annot %>%
            mutate(cell = factor(as.integer(factor(cell, cell_order)), 1:length(cell_order))) %>%
            ggplot(
                aes(x = cell, y = factor(clade), fill = I)
            ) +
            geom_tile(width=0.8, height=0.8) +
            theme_bw() +
            theme(
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = text_size),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                panel.grid = element_blank(),
                plot.margin = margin(1, 0, 0, 0, unit = "mm"), 
            ) +
            scale_x_discrete(expand = expansion(add = 1), drop = F) +
            scale_y_discrete(expand = expansion(add = 1))
    }

    if (!is.null(cell_annot)) {

        cell_annot = cell_annot %>% filter(cell %in% phylo$tip.label)

        p_bar = cell_annot %>%
            mutate(cell = factor(cell, cell_order)) %>%  
            annot_bar(legend = annot_legend, label_group = label_group, 
                label_size = text_size/2,
                annot_scale = annot_scale, legend_title = annot_legend_title)

        if (!is.null(clade_annot)) {
            (p_tree / p_bar / p_clade / p_heatmap) + plot_layout(heights = c(1, 0.1, 1, 2))
        } else {
            (p_tree / p_bar / p_heatmap) + plot_layout(heights = c(1, 0.1, 2))
        }
    } else {
        (p_tree / p_heatmap) + plot_layout(heights = c(1, 2))
    }
}


plot_phylo_circ = function(gtree, node_conf = FALSE, conf_label = FALSE, title = '', branch_width = 0.3, dot_size = 1, conf_min = 0.5, cell_annot = NULL) {

    p_tree = ggtree(gtree, layout = 'circular', branch.length = "none", linewidth = branch_width) +
            ggtitle(title)

    if (node_conf) {

        dat = gtree %>% activate(nodes) %>%
            mutate(isRoot = node_is_root()) %>%
            as.data.frame() %>% 
            select(any_of(c('name', 'isRoot', 'conf')))

        if ('conf' %in% colnames(dat)) {
            p_tree = p_tree %<+% 
                dat +
                geom_nodepoint(aes(color = conf, subset = !isTip & !isRoot, x = branch), size = dot_size, pch = 16, stroke = 1) +
                scale_color_gradient(low = 'white', high = 'firebrick', limits = c(conf_min,1), oob = scales::oob_squish)
    
            if (conf_label) {
                p_tree = p_tree + 
                    geom_text2(
                        aes(label = round(conf, 2), x = branch, subset = !isTip & !isRoot),
                        size = text_size, hjust = 0, vjust = 0.25
                    )
            }
        }
    }

    if (!is.null(cell_annot)) {
        p_tree = p_tree + geom_fruit(
                data = cell_annot,
                geom = geom_col,
                mapping = aes(y = cell, fill = annot, x = 1),
                show.legend = F
            )
    }

    return(p_tree)
}
