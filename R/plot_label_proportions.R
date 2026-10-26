#' @title Plot label proportions per group
#' @name plotLabelProportions
#' @description
#' One stacked bar per group showing how it is composed of labels: what
#' fraction of each spatial niche, cluster or region is each cell type. Groups
#' are ordered by a hierarchical clustering of their composition, with the
#' dendrogram drawn to the right of the bars.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams gmulti_params
#' @param groups character. Cell metadata column whose values form the bars,
#' such as a niche or cluster assignment.
#' @param labels character. Cell metadata column whose values are stacked
#' within each bar, such as a cell type.
#' @param spat_enr_name character. A spatial enrichment holding a composition
#' per cell, used instead of `labels`. Its columns are averaged within each
#' group, so a neighborhood composition from
#' `calculateLabelProportions(group_method = "spatialnetwork")` becomes the
#' mean neighborhood of the cells in each group.
#' @param proportions data.frame. Precomputed proportions instead of `labels`
#' or `spat_enr_name`: one row per group, the group IDs in the first column
#' and one numeric column per label, as returned by
#' `calculateLabelProportions(output = "data.table")`.
#' @param value character. "proportion" (default) scales each bar to one;
#' "count" shows the number of cells, and needs `labels`.
#' @param cluster_groups logical. Order the bars by a hierarchical clustering
#' of their composition. Otherwise they are sorted by name.
#' @param dist_method character. Distance between compositions before
#' clustering. "correlation" (1 - r) compares the shape of a composition
#' rather than its size; "euclidean" is also available.
#' @param hclust_method character. Passed to `stats::hclust()`.
#' @param show_dendrogram logical. Draw the dendrogram to the right of the
#' bars, in the same plot. Needs the ggdendro package.
#' @param custom_group_order,custom_label_order character. Explicit order of
#' the bars and of the labels within them.
#' @param label_color_code named character. A color per label, such as the
#' cell type colors used in other plots.
#' @param show_values logical. Print each segment's value on it.
#' @param min_value_shown numeric. Segments narrower than this fraction of
#' their bar are left unlabeled.
#' @details
#' Exactly one of `labels`, `spat_enr_name` and `proportions` is used. With
#' `labels` the proportions are computed by
#' `GiottoClass::calculateLabelProportions()`, so this draws what that
#' function returns. The two group answers differ: `labels` gives the share
#' of the cells in each group that carry each label, `spat_enr_name` the mean
#' composition of the neighborhoods of the cells in each group.
#' @param value_text_size numeric. Size of the values printed on the bars.
#' @returns ggplot
#' @seealso `Giotto::plotCTCompositionInNicheCluster()` for niche clusters
#' from ONTraC.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' plotLabelProportions(g, groups = "leiden_clus", labels = "custom_leiden")
#' @export
plotLabelProportions <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    groups,
    labels = NULL,
    spat_enr_name = NULL,
    proportions = NULL,
    value = c("proportion", "count"),
    cluster_groups = TRUE,
    dist_method = c("correlation", "euclidean"),
    hclust_method = "ward.D",
    show_dendrogram = TRUE,
    custom_group_order = NULL,
    custom_label_order = NULL,
    label_color_code = NULL,
    show_values = TRUE,
    min_value_shown = 0.04,
    value_text_size = 2.6,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "plotLabelProportions",
    view = NULL) {
    # NSE vars
    group <- label <- frac <- txt <- x <- y <- xend <- yend <- tx <- txend <- NULL

    value <- match.arg(value)
    dist_method <- match.arg(dist_method)
    sources <- c(labels = !is.null(labels),
        spat_enr_name = !is.null(spat_enr_name),
        proportions = !is.null(proportions))
    if (sum(sources) != 1L) {
        stop("provide exactly one of `labels`, `spat_enr_name` and ",
            "`proportions`", call. = FALSE)
    }
    if (identical(value, "count") && !sources[["labels"]]) {
        stop("value = \"count\" needs `labels`: an enrichment or a table of ",
            "proportions carries no cell counts", call. = FALSE)
    }
    checkmate::assert_flag(cluster_groups)
    checkmate::assert_flag(show_dendrogram)
    checkmate::assert_flag(show_values)
    checkmate::assert_number(min_value_shown, lower = 0, upper = 1)
    checkmate::assert_number(value_text_size, lower = 0)

    if (sources[["proportions"]]) {
        prop <- data.table::as.data.table(proportions)
        group_name <- if (missing(groups)) colnames(prop)[[1L]] else groups
    } else {
        checkmate::assert_string(groups)
        group_name <- groups
        spat_unit <- set_default_spat_unit(
            gobject = gobject, spat_unit = spat_unit)
        feat_type <- set_default_feat_type(
            gobject = gobject, spat_unit = spat_unit, feat_type = feat_type)
        gobject <- .gg_resolve(gobject, view,
            slots = c("cell_metadata",
                if (sources[["spat_enr_name"]]) "spatial_enrichment"),
            spat_unit = spat_unit, feat_type = feat_type)
        meta <- getCellMetadata(gobject, spat_unit = spat_unit,
            feat_type = feat_type, output = "data.table")
        absent <- setdiff(c(groups, labels), colnames(meta))
        if (length(absent)) {
            stop("not in the cell metadata: ", paste(absent, collapse = ", "),
                call. = FALSE)
        }
    }

    if (sources[["labels"]]) {
        prop <- calculateLabelProportions(gobject,
            labels = labels, group_method = "table", groups = groups,
            spat_unit = spat_unit, feat_type = feat_type,
            output = "data.table", verbose = FALSE)
    } else if (sources[["spat_enr_name"]]) {
        enr <- getSpatialEnrichment(gobject, spat_unit = spat_unit,
            feat_type = feat_type, name = spat_enr_name,
            output = "data.table")
        num <- setdiff(names(which(vapply(enr, is.numeric, logical(1)))),
            "cell_ID")
        if (any(enr[, num, with = FALSE] < 0, na.rm = TRUE)) {
            stop("spatial enrichment '", spat_enr_name, "' has negative ",
                "values, so it is not a composition to stack", call. = FALSE)
        }
        grp <- meta[[groups]][match(enr$cell_ID, meta$cell_ID)]
        keep <- !is.na(grp)
        prop <- enr[keep, lapply(.SD, mean, na.rm = TRUE),
            by = list(group = as.character(grp[keep])), .SDcols = num]
    }

    m <- as.matrix(prop[, -1L, with = FALSE])
    m[is.na(m)] <- 0
    rownames(m) <- as.character(prop[[1L]])
    if (identical(value, "count")) {
        sizes <- table(as.character(meta[[groups]]))
        m <- round(m * as.vector(sizes[rownames(m)]))
    } else {
        m <- m / rowSums(m)
    }

    hc <- NULL
    if (!is.null(custom_group_order)) {
        absent <- setdiff(custom_group_order, rownames(m))
        if (length(absent)) {
            stop("custom_group_order names absent groups: ",
                paste(absent, collapse = ", "), call. = FALSE)
        }
        group_order <- custom_group_order
    } else if (cluster_groups && nrow(m) > 2L && ncol(m) > 1L) {
        d <- if (identical(dist_method, "correlation")) {
            suppressWarnings(stats::as.dist(1 - stats::cor(t(m))))
        } else {
            stats::dist(m)
        }
        # a perfectly even composition has no correlation with anything
        if (anyNA(d)) d <- stats::dist(m)
        hc <- stats::hclust(d, method = hclust_method)
        group_order <- rownames(m)[hc$order]
    } else {
        group_order <- sort(rownames(m))
    }
    label_order <- if (!is.null(custom_label_order)) {
        custom_label_order
    } else {
        colnames(m)[order(-colSums(m))]
    }

    long <- data.table::as.data.table(as.table(m))
    data.table::setnames(long, c("group", "label", "value"))
    long[, `:=`(group = factor(group, levels = group_order),
        label = factor(label, levels = label_order))]

    colors <- if (!is.null(label_color_code)) {
        label_color_code
    } else {
        stats::setNames(set_default_color_discrete_cell(
            instrs = instructions(gobject))(n = length(label_order)),
            label_order)
    }

    pl <- ggplot2::ggplot(long, ggplot2::aes(value, group, fill = label)) +
        ggplot2::geom_col(width = 0.8) +
        ggplot2::scale_fill_manual(values = colors, breaks = label_order,
            na.value = "grey85") +
        ggplot2::scale_y_discrete(expand = c(0, 0.6)) +
        ggplot2::labs(y = group_name,
            fill = labels %null% spat_enr_name %null% "label",
            x = if (identical(value, "count")) "cells" else
                if (sources[["spat_enr_name"]]) "mean neighborhood composition"
                else "fraction of cells") +
        ggplot2::theme_classic() +
        ggplot2::theme(plot.margin = ggplot2::margin(5.5, 14, 5.5, 5.5))
    bar_end <- max(long[, sum(value), by = "group"]$V1)
    bar_breaks <- function(lims) {
        b <- scales::extended_breaks()(c(0, bar_end))
        b[b >= 0 & b <= bar_end]
    }
    pl <- pl + if (identical(value, "count")) {
        ggplot2::scale_x_continuous(expand = c(0, 0), breaks = bar_breaks,
            labels = function(z) format(z, big.mark = ",", trim = TRUE))
    } else {
        ggplot2::scale_x_continuous(expand = c(0, 0), breaks = bar_breaks,
            labels = function(z) paste0(z * 100, "%"))
    }

    if (show_values) {
        as_pct <- identical(value, "proportion")
        lab_dt <- data.table::copy(long)
        # counts share one axis across bars of different lengths, so a
        # segment's room is measured against the longest bar
        if (as_pct) {
            lab_dt[, frac := value / sum(value), by = "group"]
        } else {
            lab_dt[, frac := value / max(lab_dt[, sum(value), by = "group"]$V1)]
        }
        lab_dt[, txt := if (as_pct) {
            paste0(round(value * 100), "%")
        } else {
            format(value, big.mark = ",", trim = TRUE)
        }]
        # blank the text rather than drop the row: position_stack stacks
        # whatever rows it is given, so a dropped row shifts every label
        lab_dt[frac < min_value_shown * pmax(1, nchar(txt) / 3), txt := ""]
        pl <- pl + ggplot2::geom_text(data = lab_dt,
            ggplot2::aes(x = value, y = group, label = txt, group = label),
            position = ggplot2::position_stack(vjust = 0.5),
            size = value_text_size, inherit.aes = FALSE)
    }

    if (show_dendrogram && !is.null(hc)) {
        if (!requireNamespace("ggdendro", quietly = TRUE)) {
            message("ggdendro is not installed, so the dendrogram is not ",
                "drawn; the bars keep its order")
        } else {
            # the tree shares the bars' y positions, drawn to their right in
            # the same panel, so legends, labs() and theme() act on one plot
            seg <- data.table::as.data.table(ggdendro::dendro_data(
                stats::as.dendrogram(hc), type = "rectangle")$segments)
            gap <- 0.02 * bar_end
            scale_h <- 0.2 * bar_end / max(seg$y)
            seg[, `:=`(tx = bar_end + gap + y * scale_h,
                txend = bar_end + gap + yend * scale_h)]
            pl <- pl +
                ggplot2::geom_segment(data = seg,
                    ggplot2::aes(x = tx, xend = txend, y = x, yend = xend),
                    linewidth = 0.3, color = "grey35", inherit.aes = FALSE) +
                ggplot2::annotate("segment", x = 0, xend = bar_end,
                    y = -Inf, yend = -Inf) +
                ggplot2::theme(axis.line.x = ggplot2::element_blank())
        }
    }

    plot_output_handler(gobject = gobject, plot_object = pl,
        save_plot = save_plot, return_plot = return_plot,
        show_plot = show_plot, default_save_name = default_save_name,
        save_param = save_param, else_return = NULL)
}
