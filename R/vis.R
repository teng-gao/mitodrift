#' @export 
plot_phylo_heatmap2 = function(gtree, df_var, branch_width = 0.25, root_edge = TRUE, dot_size = 1, ylim = NULL,
    clade_annot = NULL, tip_annot = NULL,
    title = NULL, auc = FALSE, clone_bar = FALSE, label_site = FALSE, cell_annot = NULL, tip_lab = FALSE, node_lab = FALSE, layered = FALSE,
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
            plot.margin = margin(0, 0, 0, 0, unit = "mm"), 
            axis.title.x = element_blank(), 
            axis.line.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.x = element_blank(),
            axis.text.y = element_text(), 
            axis.line.y = element_line(),
            axis.ticks.y = element_line(), 
            axis.ticks.length.x = unit(0, "pt"),
            panel.background = element_rect(fill = "transparent", colour = NA), 
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.grid = element_blank()
        ) + 
        coord_flip() +
        scale_x_reverse(expand = expansion(mult = 0.05)) +
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
            ggraph::geom_node_point(aes(color = p_v), size = dot_size, stroke = 0) +
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

    if (!is.null(tip_annot)) {

        dat = tip_annot %>% select(any_of(c('name' = 'cell', 'annot')))

        p_tree = p_tree %<+% 
            dat +
            geom_tippoint(aes(color = annot), size = dot_size, pch = 19)
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

    if (!'vaf' %in% colnames(df_var)) {
        df_var = df_var %>% mutate(vaf = a/d)
    }

    cell_order = p_tree$data %>% filter(isTip) %>% arrange(y) %>% 
        pull(label)

    mut_order = order_muts(cell_order, df_var)
    
    df_var = df_var %>% filter(cell %in% cell_order) %>%
        mutate(variant = factor(variant, rev(mut_order))) %>%
        mutate(cell = factor(as.integer(factor(cell, cell_order)), 1:length(cell_order)))

    p_heatmap = df_var %>% 
        # filter(vaf > 0) %>%
        ggplot(
            aes(x = cell, y = variant, fill = vaf)
        ) +
        geom_raster() +
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
        guides(fill = guide_colorbar(title = 'VAF')) +
        xlab(paste0('cells (n=', length(cell_order), ')')) +
        ylab(paste0('variants (n=', length(unique(df_var$variant)), ')'))

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
                annot_scale = annot_scale, legend_title = annot_legend_title, layered = layered)

        if (!is.null(clade_annot)) {
            (p_tree / p_bar / p_clade / p_heatmap) + plot_layout(heights = c(1, 0.1, 1, 2))
        } else {
            (p_tree / p_bar / p_heatmap) + plot_layout(heights = c(1, 0.1, 2))
        }
    } else {
        (p_tree / p_heatmap) + plot_layout(heights = c(1, 2))
    }
}

# expect columns cell and annot
#' @keywords internal
annot_bar = function(
    D, transpose = FALSE, legend = TRUE, legend_title = '', size = 0.05, label_group = FALSE, label_size = 5,
    pal_annot = NULL, annot_scale = NULL, raster = FALSE, layered = FALSE
) {

    if (layered) {
        p = ggplot(D, aes(x = cell, y = annot, fill = annot)) +
            scale_y_discrete(expand = expansion(add = 1.5))
    } else {
        p = ggplot(D, aes(x = cell, y = legend_title, fill = annot)) +
            scale_y_discrete(expand = expansion(mult = 0.01))
    }

    p = p +
        geom_tile(width=0.8, height=0.8, size = 0) +
        theme_void() +
        scale_x_discrete(expand = expansion(add = 1), drop = F) +
        theme(
            panel.spacing = unit(0.1, 'mm'),
            panel.border = element_rect(size = 0, color = 'black', fill = NA),
            panel.background = element_rect(fill = 'white'),
            strip.background = element_blank(),
            strip.text = element_blank(),
            # axis.text = element_text(size = 8),
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(0.1,0,0.1,0, unit = 'mm')
        )

    if (label_group) {
        p = p + geom_text(aes(label = annot), size = label_size, angle = 90)
    }

    # if (!is.null(annot_scale)) {
    #     p = p + annot_scale
    # } else {
    #     if (is.null(pal_annot)) {
    #         pal = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    #         getPalette = colorRampPalette(pal)
    #         pal_annot = getPalette(length(unique(D$annot)))
    #     }
    #     p = p + scale_fill_manual(values = pal_annot, na.value = 'gray90', limits = force)
    # }

    p = p + scale_fill_discrete(na.value = 'gray90', limits = force)

    if (transpose) {
        p = p + coord_flip() +
            theme(plot.margin = margin(0,0.5,0,0.5, unit = 'mm'))
    }

    if (legend) {
        p = p + guides(fill = guide_legend(keywidth = unit(3, 'mm'), keyheight = unit(1, 'mm'), title = legend_title))
    } else {
        p = p + guides(fill = 'none')
    }

    if (raster) {
        p = ggrastr::rasterize(p, layers = 'Tile', dpi = 300)
    }

    return(p)
}


order_muts_binary = function(cell_order, mut_dat) {
    
    vaf_mat = mut_dat %>% reshape2::dcast(variant ~ cell, value.var = 'vaf') %>%
        tibble::column_to_rownames('variant')
    
    pres <- vaf_mat > 0
    tip_idx   <- seq_along(cell_order)
    names(tip_idx) <- cell_order
    
    # weighted-average tip for each variant
    avg_tip_pos <- sapply(rownames(vaf_mat), function(v) {
      cells <- which(pres[v,])
      mean(tip_idx[ colnames(vaf_mat)[cells] ])
    })
    
    mut_order <- names(sort(avg_tip_pos))

    return(mut_order)
}

order_muts <- function(cell_order, mut_dat) {
  
  # 1) cast to a variant × cell matrix of continuous VAFs
  vaf_mat <- mut_dat %>% 
    reshape2::dcast(variant ~ cell, value.var = 'vaf') %>%
    tibble::column_to_rownames('variant')
  
  # 2) make sure columns follow the tree’s tip order
  vaf_mat <- vaf_mat[, cell_order, drop = FALSE]
  
  # 3) create an index for each tip
  tip_idx <- seq_along(cell_order)
  names(tip_idx) <- cell_order
  
  # 4) compute weighted-median tip for each variant
  med_tip_pos <- apply(vaf_mat, 1, function(vafs) {
    # treat missing as zero
    vafs[is.na(vafs)] <- 0
    tot <- sum(vafs)
    if (tot == 0) {
      return(NA_real_)    # no signal → put at end
    }
    # cumulative weight
    cs <- cumsum(vafs)
    # first tip where cum-weight ≥ half total
    k  <- which(cs >= tot/2)[1]
    tip_idx[k]
  })
  
  # 5) sort variants by that median tip
  names(sort(med_tip_pos, na.last = TRUE))
}

plot_phylo_circ = function(gtree, node_conf = FALSE, conf_label = FALSE, title = '', pwidth = 0.25,
    branch_width = 0.3, dot_size = 1, conf_min = 0.5, cell_annot = NULL, offset = 0.05, width = 0.8,
    activity_mat = NULL,
    tip_annot = NULL, legend = FALSE, layered = FALSE) {

    p_tree = ggtree(gtree, layout = 'circular', branch.length = "none", linewidth = branch_width) +
            ggtitle(title)

    if (node_conf & inherits(gtree, 'tbl_graph')) {

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
        if (layered) {
            p_tree = p_tree + geom_fruit(
                data = cell_annot,
                geom = geom_tile,
                pwidth = pwidth,
                mapping = aes(y = cell, fill = annot, x = annot),
                show.legend = legend)
        } else {
            p_tree = p_tree + geom_fruit(
                data = cell_annot,
                geom = geom_col,
                pwidth = pwidth,
                mapping = aes(y = cell, fill = annot, x = 1),
                show.legend = legend)
        }
    }

    if (!is.null(tip_annot)) {

        dat = tip_annot %>% select(any_of(c('name' = 'cell', 'annot')))

        p_tree = p_tree %<+% 
            dat +
            geom_tippoint(aes(color = annot), size = dot_size, pch = 19, stroke = 0)
    }

    if (!is.null(activity_mat)) {

        p_tree = p_tree + geom_fruit(
            data = activity_mat %>%
                reshape2::melt() %>% 
                setNames(c('gene', 'cell', 'value')) %>%
                group_by(gene) %>%
                mutate(value = as.vector(scale(value))) %>%
                ungroup(),
            geom = geom_tile,
            offset = offset,
            width = width,
            # linewidth = 0.3,
            axis.params = list(
            axis       = "x", 
            text.angle = 45,
            text.size  = 2,
            vjust      = 0
            ),
            mapping = aes(x = gene, y = cell, fill = value),
            show.legend = TRUE
        ) +
        scale_fill_gradient2(low = "blue", high = "red", limits = c(-2,2), oob = scales::oob_squish)

        if (!is.null(tip_annot)) {
            cts = unique(tip_annot$annot)
            annot_cols = RColorBrewer::brewer.pal(length(cts), 'Set1') %>% setNames(cts)

            p_tree = p_tree +
                scale_color_manual(
                    name   = "Cell type",   
                    values = annot_cols
                ) +
                geom_point(
                    data       = data.frame(annot = names(annot_cols)),
                    aes(colour  = annot),
                    x          = 0,      # plotted off the panel
                    y          = 0,
                    size = 0, stroke = 0,
                    show.legend= TRUE
                )  +
                guides(color = guide_legend(override.aes = list('size' = 1))) 
        }
    }

    return(p_tree)
}
