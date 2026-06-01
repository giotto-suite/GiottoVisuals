#' @name dotPlot
#' @title Create a dotplot
#' @description Visualize feature expression statistics applied across
#'   clusters/groupings of cells. The default behavior is dot size scaled by
#'   what percentage of cells within a particular cluster express the feature,
#'   and dot color scaled by mean expression of that feature within the cluster.
#' @inheritParams data_access_params
#' @inheritParams plot_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_output_params
#' @param dot_size,dot_color summary function e.g. `sum`, `mean`, `var`, or
#'   other custom function. The default for `dot_size` finds the percentage of
#'   cells of a particular cluster that do not have an expression level of 0.
#' @param dot_size_threshold numeric. The minimal value at which a dot is no
#'   longer drawn.
#' @param feats character vector or named list. Features to use or named lists
#' of features to use.
#' @param cluster_column character. Clusterings column to use (usually in cell
#'   metadata)
#' @param cluster_custom_order character vector. Specific cluster order to use
#' @param dot_scale numeric. Controls size of dots
#' @param dot_color_gradient hex codes or palette name. Color gradient to use.
#' @param gradient_limits numeric vector of length 2. Set minmax value mappings
#'   for color gradient
#' @param expression_values character. Expression values to use.
#' @param title character. title for plot
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' library(GiottoClass)
#' f1 <- head(featIDs(g), 20)
#' f2 <- tail(featIDs(g), 20)
#'
#' dotPlot(g, cluster_column = "leiden_clus", feats = f1)
#'
#' # make two plots from two named sets of features
#' dotPlot(g,
#'     cluster_column = "leiden_clus",
#'     feats = list(set1 = f1, set2 = f2)
#' )
#'
#' # add new cell metadata col that classifies by total expression into
#' # - 0 (< 900)
#' # - 1 (> 900 and < 1200)
#' # - 2 (> 1200)
#' total_expr <- g$total_expr
#' g$expr_class <- findInterval(total_expr, c(900, 1200))
#'
#' # Create a dotplot while splitting the values based on the above total
#' # expression classifications.
#' dotPlot(g,
#'     cluster_column = "leiden_clus",
#'     feats = list(set1 = f1, set2 = f2),
#'     group_by = "expr_class"
#' )
#'
#' # Same as before, but focusing on classifications 0 and 2
#' dotPlot(g,
#'     cluster_column = "leiden_clus",
#'     feats = list(set1 = f1, set2 = f2),
#'     group_by = "expr_class",
#'     group_by_subset = c(0, 2)
#' )
#'
#' # example with an alternate function used for `dot_color` and a different
#' # color gradient
#' dotPlot(g,
#'     cluster_column = "leiden_clus",
#'     feats = f1,
#'     dot_size = mean,
#'     dot_color = var,
#'     dot_color_gradient = c("#EEEEFF", "#333377")
#' )
#' @export
dotPlot <- function(gobject,
    feats,
    cluster_column,
    cluster_custom_order = NULL,
    dot_size = function(x) mean(x != 0) * 100,
    dot_size_threshold = 0,
    dot_scale = 6,
    dot_color = mean,
    dot_color_gradient = NULL,
    gradient_midpoint = NULL,
    gradient_style = "sequential",
    gradient_limits = NULL,
    group_by = NULL,
    group_by_subset = NULL,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c(
        "normalized",
        "scaled",
        "custom"
    ),
    title = NULL,
    show_legend = TRUE,
    legend_text = 10,
    legend_symbol_size = 2,
    background_color = "white",
    axis_text = 10,
    axis_title = 9,
    cow_n_col = NULL,
    cow_rel_h = 1,
    cow_rel_w = 1,
    cow_align = "h",
    theme_param = list(),
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "dotPlot",
    view = NULL,
    space = NULL) {
    checkmate::assert_character(cluster_column, len = 1L)
    .gg_assert_giotto_single(gobject)
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression", "spatial_enrichment"))
    if (!is.null(gradient_limits)) {
        checkmate::assert_numeric(gradient_limits, len = 2L)
    }
    title <- title %null% ""

    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    expression_values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )

    clus <- spatValues(gobject,
        spat_unit = spat_unit, feat_type = feat_type,
        feats = cluster_column, verbose = FALSE
    )
    expr <- spatValues(gobject,
        spat_unit = spat_unit, feat_type = feat_type,
        feats = unique(unlist(feats)), # unlist to get all feats
        expression_values = expression_values,
        verbose = FALSE
    )

    # combine cluster and expression info
    ann_dt <- clus[expr, on = "cell_ID"]

    common_args <- get_args_list(keep = c(
        "cluster_column", "dot_size_threshold", "dot_size", "dot_color",
        "cluster_custom_order", "gradient_limits", "gradient_midpoint",
        "dot_color_gradient", "gradient_style", "dot_scale", "theme_param",
        "legend_text", "axis_title", "axis_text", "background_color"
    ))
    common_args$instrs <- instructions(gobject)
    # `ann_dt`, `feats`, `title` need to be updated

    # split up ann_dt into lists of ann_dt if `group_by` was given
    if (is.null(group_by)) {
        adt_list <- list(ann_dt)
    } else {
        grpby <- spatValues(gobject,
            spat_unit = spat_unit, feat_type = feat_type, feats = group_by,
            verbose = FALSE
        )
        ann_dt <- ann_dt[grpby, on = "cell_ID"]
        # consider name collisions
        if (identical(group_by, cluster_column)) {
            group_by <- paste0("i.", group_by)
        }

        unique_groups <- unique(ann_dt[[group_by]])
        # subset unique groups to those selected with `group_by_subset`
        if (!is.null(group_by_subset)) {
            not_found <- group_by_subset[!group_by_subset %in% unique_groups]
            if (length(not_found) > 0) {
                message("the following subset was not found: ", not_found)
            }
            unique_groups <- unique_groups[unique_groups %in% group_by_subset]
        }

        # subset ann_dt based on the group_by
        adt_list <- lapply(unique_groups, function(ugroup) {
            ann_dt[get(group_by) == ugroup]
        })

        # update title
        title <- paste(title, unique_groups)
    }

    # deal with multiple feature sets
    if (!is.list(feats)) {
        feats <- list(feats)
    } else {
        fnames <- names(feats)
        if (length(fnames) != length(feats)) {
            stop("If `feats` is a list, all list elements must be named")
        }

        # duplicate adt_list for however many feature sets there are
        # no need to bother with subsetting feat columns
        n_adt <- length(adt_list)
        n_fts <- length(feats)
        adt_list <- rep(adt_list, n_fts)

        # duplicate feats list to match adt list
        feats <- rep(feats, n_adt)

        # update title
        title <- paste(rep(title, n_fts), fnames)
    }


    pl <- mapply(
        function(adt_i, title_i, feats_i) {
            specific_args <- list(
                ann_dt = adt_i,
                feats = feats_i,
                title = title_i
            )
            do.call(.dplot_single, args = c(specific_args, common_args))
        },
        adt_i = adt_list, title_i = title, feats_i = feats,
        SIMPLIFY = FALSE
    )

    if (length(pl) == 1L) pl <- pl[[1L]]
    if (length(pl) > 1 && !inherits(pl, "gg")) {
        pl <- cowplot::plot_grid(
            plotlist = pl,
            ncol = set_default_cow_n_col(
                cow_n_col = cow_n_col,
                nr_plots = length(pl)
            ),
            rel_heights = cow_rel_h,
            rel_widths = cow_rel_w,
            align = cow_align
        )
    }

    plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    )
}


# internals ####

.dplot_single <- function(
    ann_dt, cluster_column, dot_size, dot_color, feats, dot_size_threshold,
    cluster_custom_order, gradient_limits, gradient_midpoint,
    dot_color_gradient, gradient_style, dot_scale, theme_param,
    legend_text, title, axis_title, axis_text, background_color, instrs) {
    # NSE vars
    cluster <- feat <- color <- size <- NULL

    dsize <- ann_dt[, lapply(.SD, dot_size), .SDcols = feats,
                    by = cluster_column]
    dcol <- ann_dt[, lapply(.SD, dot_color), .SDcols = feats,
                    by = cluster_column]

    dsize <- data.table::melt(dsize,
        id.vars = cluster_column,
        measure.vars = feats,
        value.name = "size",
        variable.name = "feat"
    )
    dcol <- data.table::melt(dcol,
        id.vars = cluster_column,
        measure.vars = feats,
        value.name = "color",
        variable.name = "feat"
    )
    plot_dt <- dsize[dcol, on = c(cluster_column, "feat")]
    data.table::setnames(plot_dt, old = cluster_column, new = "cluster")

    ## dot size cutoff ##
    plot_dt <- plot_dt[size > dot_size_threshold, ]

    ## set cluster order ##
    if (is.null(cluster_custom_order)) {
        plot_dt[
            , cluster := factor(cluster, levels = mixedsort(unique(cluster)))]
    } else {
        plot_dt[, cluster := factor(cluster, levels = cluster_custom_order)]
    }

    # apply limits
    if (!is.null(gradient_limits)) {
        plot_dt[, color := scales::oob_squish(color, gradient_limits)]
    }

    pl <- ggplot2::ggplot() +
        ggplot2::geom_point(
            data = plot_dt,
            ggplot2::aes(x = cluster, y = feat, color = color, size = size)
        )

    # apply color gradient
    if (is.null(gradient_midpoint)) {
        gradient_midpoint <-
            stats::median(plot_dt$color)
    }
    pl <- pl + set_default_color_continuous_cell(
        colors = dot_color_gradient,
        instrs = instrs,
        midpoint = gradient_midpoint,
        style = gradient_style,
        type = "color"
    )

    # size scaling
    pl <- pl + scale_size_continuous(range = c(1, dot_scale))

    ## theme ##
    gg_theme_args <- c(
        theme_param,
        legend_text = legend_text,
        axis_title = axis_title,
        axis_text = axis_text,
        background_color = background_color,
        axis_text_y_angle = 0,
        axis.ticks = element_blank()
    )
    pl <- pl + do.call(.gg_theme, args = gg_theme_args)

    pl <- pl + ggplot2::labs(title = title)

    return(pl)
}
