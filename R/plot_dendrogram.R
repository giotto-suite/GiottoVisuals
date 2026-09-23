#' @title showClusterDendrogram
#' @name showClusterDendrogram
#' @description Creates dendrogram for selected clusters.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values expression values to use
#' (e.g. "normalized", "scaled", "custom")
#' @param cluster_column name of column to use for clusters
#' (e.g. "leiden_clus")
#' @param cor correlation score to calculate distance
#' (e.g. "pearson", "spearman")
#' @param distance distance method to use for hierarchical clustering,
#' default to "ward.D"
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param rotate rotate dendrogram 90 degrees
#' @param tree optional `hclust` from `Giotto::calculateClusterTree()`. When
#' supplied the tree is plotted as given rather than rebuilt, so the same tree
#' can back the plot, the splits and any per-node analysis.
#' @inheritParams gmulti_params
#' @inheritDotParams ggdendro::ggdendrogram
#' @details Expression correlation dendrogram for selected clusters.
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' showClusterDendrogram(g, cluster_column = "leiden_clus")
#'
#' @export
showClusterDendrogram <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    cluster_column,
    cor = c("pearson", "spearman"),
    distance = "ward.D",
    h = NULL,
    h_color = "red",
    rotate = FALSE,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "showClusterDendrogram",
    tree = NULL,
    view = NULL,
    space = NULL,
    ...) {
    # verify if optional package is installed
    package_check(pkg_name = "ggdendro", repository = "CRAN")
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression", "spatial_enrichment"))

    # A tree from `Giotto::calculateClusterTree()` can back the dendrogram, the
    # splits and any per-node analysis at once, instead of each rebuilding its
    # own from the expression values and being free to disagree. It also owns
    # the cluster ordering the correlation depends on.
    #
    # Placed after `.gg_materialize()`, not before it: `gobject` still reaches
    # `plot_output_handler()` below, so `view` / `space` and a `giottoMulti`
    # behave the same whether or not a tree was supplied.
    if (!is.null(tree)) {
        if (!inherits(tree, "hclust")) {
            stop("[showClusterDendrogram] `tree` must be an `hclust`, as ",
                "returned by `Giotto::calculateClusterTree()`. Got: ",
                paste(class(tree), collapse = "/"), ".", call. = FALSE)
        }
        pl <- ggdendro::ggdendrogram(
            data = stats::as.dendrogram(tree), rotate = rotate, ...
        )
        if (!is.null(h)) {
            pl <- pl + ggplot2::geom_hline(yintercept = h, col = h_color)
        }
        return(plot_output_handler(
            gobject = gobject, plot_object = pl, save_plot = save_plot,
            return_plot = return_plot, show_plot = show_plot,
            default_save_name = default_save_name, save_param = save_param,
            else_return = NULL
        ))
    }

    values <- match.arg(
        expression_values,
        unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
    )

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    metatable <- calculateMetaTable(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        metadata_cols = cluster_column
    )

    pl <- create_cluster_dendrogram(
        data = metatable,
        clus_col = "uniq_ID",
        var_col = "variable",
        val_col = "value",
        cor = cor,
        distance = distance,
        h = h,
        h_color = h_color,
        rotate = rotate,
        ...
    )

    return(plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}


#' @name create_cluster_dendrogram
#' @title Create clustered expression dendrogram
#' @description Create a dendrogram based on a data.table with columns for
#' cluster ID, variables, and their values. If no specific values are provided
#' for the 'col' params then they will be assumed as 1. clus_col, 2. var_col,
#' 3. val_col
#' @param data data.table. Should include columns with clustering labels,
#' variables that are being clustered, and the values of those clusters.
#' @param clus_col character. name of column with clustering info
#' @param var_col character. name of column with variable name
#' @param val_col character. name of column with values info
#' @param cor correlation score to calculate distance
#' (e.g. "pearson", "spearman")
#' @param distance distance method to use for hierarchical clustering,
#' default to "ward.D"
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param rotate rotate dendrogram 90 degrees
#' @inheritDotParams ggdendro::ggdendrogram
#' @returns ggdendrogram
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#'
#' g_expression <- head(GiottoClass::getExpression(g, output = "matrix"))
#' g_expression_df <- as.data.frame(as.matrix(g_expression))
#' g_expression_df$feat_ID <- rownames(g_expression)
#'
#' g_expression_melt <- data.table::melt(g_expression_df,
#'     id.vars = "feat_ID",
#'     measure.vars = colnames(g_expression), variable.name = "cell_ID",
#'     value.name = "raw_expression"
#' )
#'
#' create_cluster_dendrogram(data.table::as.data.table(g_expression_melt),
#'     var_col = "cell_ID", clus_col = "feat_ID", "raw_expression"
#' )
#'
#' @export
create_cluster_dendrogram <- function(
        data,
        clus_col = names(data)[[1]],
        var_col = names(data)[[2]],
        val_col = names(data)[[3]],
        cor = c("pearson", "spearman"),
        distance = "ward.D",
        h = NULL,
        h_color = "red",
        rotate = FALSE,
        ...) {
    checkmate::assert_data_table(data)
    checkmate::assert_character(clus_col)
    checkmate::assert_character(var_col)
    checkmate::assert_character(val_col)
    cor <- match.arg(cor, c("pearson", "spearman"))

    dcast_metatable <- data.table::dcast.data.table(
        data = data,
        formula = paste0(var_col, "~", clus_col),
        value.var = val_col
    )
    testmatrix <- dt_to_matrix(x = dcast_metatable)

    # correlation
    cormatrix <- cor_flex(x = testmatrix, method = cor)
    cordist <- stats::as.dist(1 - cormatrix, diag = TRUE, upper = TRUE)
    corclus <- stats::hclust(d = cordist, method = distance)

    cordend <- stats::as.dendrogram(object = corclus)

    # plot dendrogram
    pl <- ggdendro::ggdendrogram(data = cordend, rotate = rotate, ...)

    # add horizontal or vertical lines
    if (!is.null(h)) {
        pl <- pl + ggplot2::geom_hline(yintercept = h, col = h_color)
    }

    pl
}


#' @title plotClusterTree
#' @name plotClusterTree
#' @description Publication figure for an annotated cluster tree: the
#' dendrogram, a ladder showing how labels merge as the tree is cut more
#' coarsely, and per-cluster evidence tracks.
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param cluster_column name of the cell metadata column holding the clusters
#' @param tree an `hclust` over the clusters, from
#' `Giotto::calculateClusterTree()`
#' @param labels the annotation, as a list with `clusters` and optionally
#' `nodes`; the same object `Giotto::annotateClusterTree()` takes
#' @param k granularities to draw as ladder bands, passed to [stats::cutree()].
#' The finest is drawn at the top, nearest the leaves.
#' @param gini_markers gini marker table carrying `detection_margin`, for the
#' specificity track. Omit to leave that track out.
#' @param margin_cut detection margin below which a leaf is marked as having no
#' feature of its own
#' @param node_labels draw the label of each internal node on the tree
#' @param palette colours for the coarsest ladder band, recycled across the
#' others. A named vector keys colours to labels.
#' @param label_size,base_size text sizes
#' @param show_plot,return_plot,save_plot,save_param,default_save_name
#' plot output settings, see [plot_output_handler()]
#' @returns ggplot
#' @details
#' The figure answers a question a bare dendrogram cannot: **at what level
#' should these clusters be annotated?** The ladder makes the trade-off legible
#' — read down to see distinct labels collapse into one, read up to see a label
#' split into types that the data may or may not support — and the tracks put
#' the evidence for each leaf next to the name it was given, so a small cluster
#' with no marker of its own is visible rather than implied.
#' @seealso `Giotto::calculateClusterTree()`,
#' `Giotto::annotateClusterTree()`, [showClusterDendrogram()]
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' tree <- Giotto::calculateClusterTree(g, cluster_column = "leiden_clus")
#' labs <- list(clusters = stats::setNames(
#'     paste("type", tree$labels), tree$labels
#' ))
#' plotClusterTree(g,
#'     cluster_column = "leiden_clus", tree = tree,
#'     labels = labs, k = c(2, 4)
#' )
#' @export
plotClusterTree <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    cluster_column,
    tree,
    labels,
    k = NULL,
    gini_markers = NULL,
    margin_cut = 25,
    node_labels = TRUE,
    palette = NULL,
    label_size = 3,
    base_size = 11,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "plotClusterTree") {
    package_check(pkg_name = "ggdendro", repository = "CRAN")

    # data.table variables
    x <- NULL

    if (!inherits(tree, "hclust")) {
        stop("[plotClusterTree] `tree` must be an `hclust`, as returned by ",
            "`Giotto::calculateClusterTree()`. Got: ",
            paste(class(tree), collapse = "/"), ".", call. = FALSE)
    }
    # `cutree()` refuses an out-of-range `k`, but in its own terms rather than
    # naming the argument the caller passed. Same guard as
    # `Giotto::annotateClusterTree()`, so the two agree on what is askable.
    if (!is.null(k)) {
        if (!is.numeric(k) || anyNA(k)) {
            stop("[plotClusterTree] `k` must be numeric and not NA.",
                call. = FALSE)
        }
        bad <- k[k < 1 | k > length(tree$labels)]
        if (length(bad)) {
            stop("[plotClusterTree] `k` must be between 1 and the number of ",
                "clusters in the tree (", length(tree$labels), "). Got: ",
                paste(unique(bad), collapse = ", "), ".", call. = FALSE)
        }
    }
    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    ord <- tree$labels[tree$order]
    n_leaf <- length(ord)
    leaf_lab <- .pct_label_vector(labels$clusters)
    node_lab <- if (!is.null(labels$nodes)) {
        .pct_label_vector(labels$nodes)
    } else {
        stats::setNames(character(0L), character(0L))
    }

    # ---- per-leaf evidence ------------------------------------------------
    cm <- GiottoClass::getCellMetadata(gobject,
        spat_unit = spat_unit, feat_type = feat_type,
        output = "data.table", copy_obj = TRUE
    )
    clus <- as.character(cm[[cluster_column]])
    size <- vapply(ord, function(x) sum(clus == x), integer(1L))

    margin <- NULL
    if (!is.null(gini_markers)) {
        gm <- data.table::as.data.table(gini_markers)
        if ("detection_margin" %in% names(gm)) {
            agg <- gm[, list(m = max(.SD[["detection_margin"]], na.rm = TRUE)),
                by = "cluster"]
            margin <- stats::setNames(agg$m, as.character(agg$cluster))[ord]
        }
    }

    # ---- ladder -----------------------------------------------------------
    # Finest first, so the band nearest the leaves is the one that matches
    # them. `cutree` returns groups keyed by leaf, which is what the band
    # geometry needs -- no extra bookkeeping.
    k <- sort(unique(k), decreasing = TRUE)
    band <- NULL
    if (length(k)) {
        band <- data.table::rbindlist(lapply(seq_along(k), function(i) {
            grp <- stats::cutree(tree, k = k[i])
            lab <- .pct_group_labels(tree, leaf_lab, node_lab, grp)
            data.table::data.table(
                x = seq_len(n_leaf),
                row = i,
                level = sprintf("k = %d", k[i]),
                label = unname(lab[ord])
            )
        }))
        band[, "level" := factor(.SD[["level"]],
            levels = sprintf("k = %d", k))]
    }

    # Colour by the coarsest band, so related leaves share a hue throughout.
    #
    # Deliberately NOT named `key`: `[.data.table` evaluates `setkeyv(ans, key)`
    # internally when `by=` is given, and a local `key` in this frame is what
    # that resolves to -- the grouped call below then fails trying to key the
    # result by a vector of cell type names.
    grp_key <- if (!is.null(band)) {
        band[band$row == length(k)][order(x)][["label"]]
    } else {
        unname(leaf_lab[ord])
    }
    if (!is.null(band)) band[, "coarse" := grp_key[.SD[["x"]]]]
    ukey <- unique(grp_key)
    cols <- if (!is.null(palette)) {
        palette
    } else {
        stats::setNames(getDistinctColors(length(ukey)), ukey)
    }

    # ---- panel 1: the tree ------------------------------------------------
    dd <- ggdendro::dendro_data(tree, type = "rectangle")
    seg <- ggdendro::segment(dd)
    p_tree <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = seg,
            ggplot2::aes(x = .data$x, y = .data$y,
                xend = .data$xend, yend = .data$yend),
            colour = "grey40", linewidth = 0.4
        ) +
        ggplot2::scale_x_continuous(
            limits = c(0.5, n_leaf + 0.5), expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(
            mult = c(0.02, 0.05)
        )) +
        ggplot2::labs(x = NULL, y = "1 - correlation") +
        ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(6, 6, 0, 6)
        )

    if (isTRUE(node_labels) && length(node_lab)) {
        nd <- .pct_node_positions(tree, node_lab)
        if (nrow(nd)) {
            p_tree <- p_tree + ggplot2::geom_label(
                data = nd,
                ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                size = label_size * 0.8, colour = "grey20",
                fill = "white", linewidth = 0,
                label.padding = ggplot2::unit(0.1, "lines")
            )
        }
    }

    # ---- panel 2: the ladder ----------------------------------------------
    pieces <- list(p_tree)
    heights <- c(3)
    if (!is.null(band)) {
        # one rectangle per run of identical labels, so a group reads as a
        # single block with one centred name rather than n repeated tiles
        runs <- band[, {
            r <- rle(.SD[["label"]])
            e <- cumsum(r$lengths)
            s <- e - r$lengths + 1L
            list(xmin = s - 0.5, xmax = e + 0.5,
                mid = (s + e) / 2, label = r$values,
                coarse = .SD[["coarse"]][s])
        }, by = "level"]

        p_band <- ggplot2::ggplot(runs) +
            ggplot2::geom_rect(
                ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                    ymin = 0.1, ymax = 0.9, fill = .data$coarse),
                colour = "white", linewidth = 0.6, alpha = 0.85
            ) +
            ggplot2::geom_text(
                ggplot2::aes(x = .data$mid, y = 0.5, label = .data$label),
                size = label_size * 0.78, colour = "grey15", check_overlap = TRUE
            ) +
            ggplot2::facet_grid(rows = ggplot2::vars(.data$level),
                switch = "y") +
            ggplot2::scale_fill_manual(values = cols, guide = "none") +
            ggplot2::scale_x_continuous(
                limits = c(0.5, n_leaf + 0.5), expand = c(0, 0)
            ) +
            ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            ggplot2::labs(x = NULL, y = NULL) +
            ggplot2::theme_minimal(base_size = base_size) +
            ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                strip.text.y.left = ggplot2::element_text(
                    angle = 0, hjust = 1, size = base_size * 0.75
                ),
                panel.spacing.y = ggplot2::unit(1.5, "pt"),
                plot.margin = ggplot2::margin(0, 6, 0, 6)
            )
        pieces <- c(pieces, list(p_band))
        heights <- c(heights, 0.42 * length(k))
    }

    # ---- panel 3: leaf labels + tracks ------------------------------------
    leaf_dt <- data.table::data.table(
        x = seq_len(n_leaf),
        cluster = ord,
        label = unname(leaf_lab[ord]),
        size = as.numeric(size),
        grp_key = grp_key
    )
    leaf_dt[, "txt" := sprintf("%s  %s", .SD[["cluster"]], .SD[["label"]])]

    p_leaf <- ggplot2::ggplot(leaf_dt) +
        ggplot2::geom_text(
            ggplot2::aes(x = .data$x, y = 0, label = .data$txt,
                colour = .data$grp_key),
            angle = 90, hjust = 1, vjust = 0.5, size = label_size,
            fontface = "bold"
        ) +
        ggplot2::scale_colour_manual(values = cols, guide = "none") +
        ggplot2::scale_x_continuous(
            limits = c(0.5, n_leaf + 0.5), expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(limits = c(-1.05, 0.05),
            expand = c(0, 0)) +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = base_size) +
        ggplot2::theme(plot.margin = ggplot2::margin(2, 6, 0, 6))
    pieces <- c(pieces, list(p_leaf))
    # scaled to the longest label: a fixed height either clips long cell type
    # names or leaves a band of dead space above the tracks
    heights <- c(heights, max(1.4, 0.052 * max(nchar(leaf_dt[["txt"]]))))

    p_size <- .pct_track(leaf_dt, "size", n_leaf, "cells\n(log10)", base_size,
        log10 = TRUE)
    pieces <- c(pieces, list(p_size))
    heights <- c(heights, 0.8)

    if (!is.null(margin)) {
        leaf_dt[, "margin" := as.numeric(margin)]
        p_margin <- .pct_track(leaf_dt, "margin", n_leaf,
            "detection\nmargin", base_size, cut = margin_cut)
        pieces <- c(pieces, list(p_margin))
        heights <- c(heights, 0.8)
    }

    # cowplot rather than patchwork: it is already an import, and `align`
    # with `axis` is what makes the panels share a leaf position -- every
    # panel is drawn on the same continuous x, so a mis-alignment here would
    # silently put a label under the wrong leaf.
    pl <- cowplot::plot_grid(
        plotlist = pieces, ncol = 1, align = "v", axis = "lr",
        rel_heights = heights
    )

    plot_output_handler(
        gobject = gobject, plot_object = pl, save_plot = save_plot,
        return_plot = return_plot, show_plot = show_plot,
        default_save_name = default_save_name, save_param = save_param,
        else_return = NULL
    )
}


# One evidence track. Bars rather than a heat strip: a length is read
# quantitatively where a colour is not, and these are the numbers a reader is
# meant to judge a label against.
#' @keywords internal
#' @noRd
.pct_track <- function(dt, col, n_leaf, lab, base_size, log10 = FALSE,
    cut = NULL) {
    d <- data.table::copy(dt)
    d[, "v" := as.numeric(.SD[[col]])]
    if (isTRUE(log10)) d[, "v" := log10(pmax(.SD[["v"]], 1))]
    d[, "flag" := if (is.null(cut)) FALSE else .SD[["v"]] < cut]

    p <- ggplot2::ggplot(d) +
        ggplot2::geom_col(
            ggplot2::aes(x = .data$x, y = .data$v, fill = .data$flag),
            width = 0.75
        ) +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "grey45", "TRUE" = "#C4622D"), guide = "none"
        ) +
        ggplot2::scale_x_continuous(
            limits = c(0.5, n_leaf + 0.5), expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(
            mult = c(0, 0.05)
        ), n.breaks = 3) +
        ggplot2::labs(x = NULL, y = lab) +
        ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = base_size * 0.6),
            axis.title.y = ggplot2::element_text(size = base_size * 0.7),
            plot.margin = ggplot2::margin(0, 6, 0, 6)
        )
    if (!is.null(cut)) {
        p <- p + ggplot2::geom_hline(yintercept = cut, linetype = 2,
            colour = "grey40", linewidth = 0.3)
    }
    p
}


# Accept a named vector or a two-column table, so `fromJSON()` output can be
# handed over unreshaped. Mirrors `Giotto:::.tree_label_vector()`, but this
# package cannot reach into Giotto's internals.
#' @keywords internal
#' @noRd
.pct_label_vector <- function(x) {
    if (is.null(x)) return(stats::setNames(character(0L), character(0L)))
    if (is.data.frame(x)) {
        key <- intersect(c("cluster", "node", "id"), names(x))[1L]
        val <- intersect(c("cell_type", "label", "name"), names(x))[1L]
        if (is.na(key) || is.na(val)) {
            stop("[plotClusterTree] label tables need an id column ",
                "(cluster/node/id) and a label column (cell_type/label/name). ",
                "Got: ", paste(names(x), collapse = ", "), ".", call. = FALSE)
        }
        x <- stats::setNames(as.character(x[[val]]), as.character(x[[key]]))
    }
    stats::setNames(as.character(x), as.character(names(x)))
}


# Label per cut group, by the same resolution order
# `Giotto::annotateClusterTree()` documents: own node, nearest labelled
# ancestor, then the leaves themselves. Kept in step with that function --
# a ladder that disagreed with the column written onto the object would be
# worse than no ladder.
#' @keywords internal
#' @noRd
.pct_group_labels <- function(tree, leaf_lab, node_lab, grp) {
    leaves_of <- function(node) {
        if (node < 0L) return(tree$labels[-node])
        c(leaves_of(tree$merge[node, 1L]), leaves_of(tree$merge[node, 2L]))
    }
    sets <- lapply(seq_len(nrow(tree$merge)), function(i) sort(leaves_of(i)))
    keys <- vapply(sets, paste, character(1L), collapse = "\r")

    out <- character(0L)
    for (gi in unique(grp)) {
        members <- sort(names(grp)[grp == gi])
        lab <- NA_character_
        hit <- match(paste(members, collapse = "\r"), keys)
        if (!is.na(hit) && as.character(hit) %in% names(node_lab)) {
            lab <- node_lab[[as.character(hit)]]
        }
        # a singleton is its own leaf, checked before any ancestor: no
        # internal node spans one leaf, so the finest band would otherwise
        # inherit a clade name and disagree with the leaf labels below it
        if (is.na(lab) && length(members) == 1L && members %in% names(leaf_lab)) {
            lab <- leaf_lab[[members]]
        }
        if (is.na(lab) && length(node_lab)) {
            contains <- vapply(sets, function(s) all(members %in% s),
                logical(1L))
            cand <- which(contains &
                as.character(seq_along(sets)) %in% names(node_lab))
            if (length(cand)) {
                lab <- node_lab[[as.character(
                    cand[which.min(lengths(sets)[cand])]
                )]]
            }
        }
        if (is.na(lab)) {
            lv <- leaf_lab[members]
            lv <- lv[!is.na(lv)]
            lab <- if (length(lv)) {
                names(sort(table(lv), decreasing = TRUE))[1L]
            } else {
                paste(members, collapse = ",")
            }
        }
        out[members] <- lab
    }
    out
}


# x/y for each labelled internal node: x is the mean position of its children,
# y its merge height, matching where `ggdendro` draws the join.
#' @keywords internal
#' @noRd
.pct_node_positions <- function(tree, node_lab) {
    pos <- match(tree$labels, tree$labels[tree$order])
    node_x <- function(node) {
        if (node < 0L) return(pos[-node])
        mean(c(node_x(tree$merge[node, 1L]), node_x(tree$merge[node, 2L])))
    }
    ids <- intersect(as.character(seq_len(nrow(tree$merge))), names(node_lab))
    if (!length(ids)) {
        return(data.table::data.table(
            x = numeric(0L), y = numeric(0L), label = character(0L)
        ))
    }
    i <- as.integer(ids)
    data.table::data.table(
        x = vapply(i, node_x, numeric(1L)),
        y = tree$height[i],
        label = unname(node_lab[ids])
    )
}
