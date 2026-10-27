#' @title Plot label composition per group
#' @name plotLabelComposition
#' @description
#' One stacked bar per group showing how it is composed of labels: what
#' fraction of each spatial niche, cluster or region is each cell type. Groups
#' are ordered by a hierarchical clustering of their composition, with the
#' dendrogram drawn to the right of the bars.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams gmulti_params
#' @param gobject giotto object. May be NULL with `proportions`; without one,
#' the plot cannot be saved and the default colors come from
#' `giotto.color_d_pal`.
#' @param groups character. Cell metadata column whose values form the bars,
#' such as a niche or cluster assignment. With `proportions`, only the axis
#' title (default: the table's first column name).
#' @param labels character. One of three inputs (see Details). Cell metadata
#' column whose values are stacked within each bar, such as a cell type.
#' @param spat_enr_name character. One of three inputs (see Details). A
#' spatial enrichment holding a composition per cell. Its columns are averaged
#' within each group, so a neighborhood composition from
#' `calculateLabelProportions(group_method = "spatialnetwork")` becomes the
#' mean neighborhood of the cells in each group. Values must be non-negative,
#' each numeric column becomes one stacked label, and `value = "count"` is not
#' available with it.
#' @param proportions data.frame. One of three inputs (see Details).
#' Precomputed proportions: one row per group, the group IDs in the first
#' column and one numeric column per label, as returned by
#' `calculateLabelProportions(output = "data.table")`. Each row is rescaled to
#' sum to one, and missing values count as zero.
#' @param value character. "proportion" (default) draws each bar to 100%;
#' "count" shows the number of cells, and needs `labels`.
#' @param cluster_groups logical. Order the bars by a hierarchical clustering
#' of their composition. Otherwise they are sorted by name. Ignored when
#' `custom_group_order` is given. With fewer than three groups or one label,
#' bars are sorted by name.
#' @param dist_method character. Distance between compositions before
#' clustering. "correlation" (1 - r) compares the shape of a composition
#' rather than its size; "euclidean" is also available. Falls back to
#' Euclidean when a correlation is undefined, as for a group with an even
#' composition.
#' @param hclust_method character. Passed to `stats::hclust()`.
#' @param show_dendrogram logical. Draw the dendrogram to the right of the
#' bars, in the same plot. Only when the bars are clustered. Without
#' ggdendro, the order is kept and the tree is not drawn.
#' @param custom_group_order,custom_label_order character. Explicit order of
#' the bars and of the labels within them. Must list every group or label in
#' the data, and nothing else.
#' @param label_color_code named character. A color per label, such as the
#' cell type colors used in other plots. Labels not named here are drawn
#' gray. Default: the suite discrete palette (`giotto.color_d_pal`).
#' @param show_values logical. Print each segment's percentage (or cell count)
#' inside it.
#' @param min_value_shown numeric. Smallest segment that gets a printed value,
#' as a fraction of its bar (of the longest bar when `value = "count"`).
#' Longer labels need proportionally more room.
#' @param value_text_size numeric. Size of the values printed on the bars.
#' @details
#' Exactly one of `labels`, `spat_enr_name` and `proportions` is used. With
#' `labels` the proportions are computed by
#' `GiottoClass::calculateLabelProportions()`, so this draws what that
#' function returns. The two group answers differ: `labels` gives the share
#' of the cells in each group that carry each label, `spat_enr_name` the mean
#' composition of the neighborhoods of the cells in each group.
#' @returns ggplot
#' @seealso `Giotto::plotCTCompositionInNicheCluster()` for niche clusters
#' from ONTraC.
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' plotLabelComposition(g, groups = "leiden_clus", labels = "custom_leiden")
#' @export
plotLabelComposition <- function(gobject = NULL,
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
    default_save_name = "plotLabelComposition",
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
    if (is.null(gobject) && !sources[["proportions"]]) {
        stop("`gobject` is needed unless `proportions` is given",
            call. = FALSE)
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

    .lp_check_order(custom_group_order, rownames(m),
        "custom_group_order", "groups")
    .lp_check_order(custom_label_order, colnames(m),
        "custom_label_order", "labels")

    hc <- NULL
    if (!is.null(custom_group_order)) {
        group_order <- custom_group_order
    } else if (cluster_groups && nrow(m) > 2L && ncol(m) > 1L) {
        d <- if (identical(dist_method, "correlation")) {
            suppressWarnings(stats::as.dist(1 - stats::cor(t(m))))
        } else {
            stats::dist(m)
        }
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
            instrs = if (!is.null(gobject)) instructions(gobject)
        )(n = length(label_order)),
            label_order)
    }

    pl <- ggplot2::ggplot(long, ggplot2::aes(value, group, fill = label)) +
        ggplot2::geom_col(width = 0.8) +
        ggplot2::scale_fill_manual(values = colors, breaks = label_order,
            na.value = "gray85") +
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
        # blank rather than drop: position_stack would shift the other labels
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
            seg <- data.table::as.data.table(ggdendro::dendro_data(
                stats::as.dendrogram(hc), type = "rectangle")$segments)
            gap <- 0.02 * bar_end
            scale_h <- 0.2 * bar_end / max(seg$y)
            seg[, `:=`(tx = bar_end + gap + y * scale_h,
                txend = bar_end + gap + yend * scale_h)]
            pl <- pl +
                ggplot2::geom_segment(data = seg,
                    ggplot2::aes(x = tx, xend = txend, y = x, yend = xend),
                    linewidth = 0.3, color = "gray35", inherit.aes = FALSE) +
                ggplot2::annotate("segment", x = 0, xend = bar_end,
                    y = -Inf, yend = -Inf) +
                ggplot2::theme(axis.line.x = ggplot2::element_blank())
        }
    }

    # plot_output_handler() needs a gobject for its show/save defaults
    if (is.null(gobject)) {
        if (isTRUE(save_plot)) {
            warning("save_plot needs a gobject to read its save settings ",
                "from", call. = FALSE)
        }
        if (!isFALSE(show_plot)) print(pl)
        return(if (isFALSE(return_plot)) invisible(NULL) else pl)
    }

    plot_output_handler(gobject = gobject, plot_object = pl,
        save_plot = save_plot, return_plot = return_plot,
        show_plot = show_plot, default_save_name = default_save_name,
        save_param = save_param, else_return = NULL)
}

.lp_check_order <- function(custom_order, present, arg, what) {
    if (is.null(custom_order)) return(invisible())
    unknown <- setdiff(custom_order, present)
    if (length(unknown)) {
        stop(arg, " names ", what, " not in the data: ",
            paste(unknown, collapse = ", "), call. = FALSE)
    }
    left_out <- setdiff(present, custom_order)
    if (length(left_out)) {
        stop(arg, " is missing ", what, ": ",
            paste(left_out, collapse = ", "), call. = FALSE)
    }
}
