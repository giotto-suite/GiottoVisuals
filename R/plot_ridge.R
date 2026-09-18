#' @name ridgePlot
#' @title Create a ridgeline plot
#' @description
#' Plot continuous numeric features as a series of density ridges. Particularly
#' useful for making comparisons between feature distributions within a
#' population.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_params
#' @inheritParams plot_cow_params
#' @param feats features to plot. Passed to [GiottoClass::spatValues()]
#' @param ... additional params to pass to [GiottoClass::spatValues()]
#' @param stat character. A `stat` accepted by
#' [ggridges::geom_density_ridges()]
#' @param ridge_alpha numeric. Alpha to use for ridgeline plot. Value from 0
#' to 1.
#' @param ridge_scale numeric. Controls height of ridges.
#' @param scale_axis character (default = `"log2"`.
#' An x axis scaling to use. This value is passed to
#' [ggplot2::scale_x_continuous()]. Set `"identity"` to remove scaling.
#' @param axis_offset numeric. Value to add to values. When this value is not
#' provided and `axis_scale` is one of `"log"`, `"log2"`, or `"log10"`, a
#' `axis_offset` of 1 is used by default.
#' @param background_color background color
#' @param show_legend logical. Whether to show legend
#' @param remove_na logical. If `FALSE`, the default, missing values are
#' removed with a warning. If `TRUE`, missing values are silently removed.
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @param title character. Title for plot
#' @param verbose be verbose
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' # plot normalized expresion of first 6 features
#' ridgePlot(g,
#'     feats = head(rownames(g)),
#'     expression_values = "normalized",
#'     scale_axis = "log2"
#' )
#'
#' # plot distribution on PCA PCs
#' ridgePlot(g,
#'     feats = sprintf("Dim.%d", 1:20),
#'     scale_axis = "identity"
#' )
#'
#' ridgePlot(g,
#'     feats = sprintf("Dim.%d", 1:20),
#'     scale_axis = "identity",
#'     group_by = "leiden_clus",
#'     group_by_subset = c(1:4)
#' )
#' @export
ridgePlot <- function(gobject, feats, ...,
    group_by = NULL,
    group_by_subset = NULL,
    spat_unit = NULL,
    feat_type = NULL,
    stat = "density_ridges",
    ridge_alpha = 1,
    ridge_scale = 5,
    scale_axis = "identity",
    axis_offset = NULL,
    background_color = "white",
    remove_na = NULL,
    show_legend = FALSE,
    title = NULL,
    theme_param = list(),
    cow_n_col = NULL,
    cow_rel_h = 1,
    cow_rel_w = 1,
    cow_align = "h",
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "ridgePlot",
    verbose = NULL,
    view = NULL,
    space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression",
            "spatial_enrichment", "dimension_reduction"))
    value <- NULL # NSE vars
    package_check("ggridges")
    checkmate::assert_character(feats)
    checkmate::assert_character(scale_axis)
    checkmate::assert_numeric(axis_offset, null.ok = TRUE)
    if (is.null(axis_offset) && scale_axis %in% c("log", "log2", "log10")) {
        axis_offset <- 1
    }
    title <- title %null% ""

    spat_unit <- set_default_spat_unit(
        gobject = gobject, spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject, spat_unit = spat_unit, feat_type = feat_type
    )

    # pull data to plot
    data <- spatValues(gobject,
        feats = feats,
        spat_unit = spat_unit,
        feat_type = feat_type,
        ...,
        verbose = verbose
    )

    # split up data into lists of data if `group_by` was given
    if (is.null(group_by)) {
        data_list <- list(
            .ridgeplot_data_prep(data, axis_offset = axis_offset)
        )
    } else {
        grpby <- spatValues(gobject,
            spat_unit = spat_unit, feat_type = feat_type, feats = group_by,
            verbose = FALSE
        )
        data <- data[grpby, on = "cell_ID"]
        unique_groups <- unique(data[[group_by]])
        # subset unique groups to those selected with `group_by_subset`
        if (!is.null(group_by_subset)) {
            not_found <- group_by_subset[!group_by_subset %in% unique_groups]
            if (length(not_found) > 0) {
                message("the following subset was not found: ", not_found)
            }
            unique_groups <- unique_groups[unique_groups %in% group_by_subset]
        }

        # subset data based on the group_by
        data_list <- lapply(unique_groups, function(ugroup) {
            data <- data[get(group_by) == ugroup]
            data <- data[, -group_by, with = FALSE]
            .ridgeplot_data_prep(data, axis_offset = axis_offset)
        })

        # update title
        title <- paste(title, unique_groups)
    }

    common_args <- get_args_list(keep = c(
        "scale_axis", "ridge_alpha", "ridge_scale", "stat", "remove_na",
        "background_color", "show_legend", "theme_param"
    ))

    pl <- mapply(
        function(data_i, title_i) {
            specific_args <- list(
                data = data_i,
                title = title_i
            )
            do.call(.ridgeplot_single, args = c(specific_args, common_args))
        },
        data_i = data_list, title_i = title,
        SIMPLIFY = FALSE
    )

    if (length(pl) == 1L) pl <- pl[[1L]] # unlist if only length 1
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

# transform table to long format and add any offset. drops cell_ID col
.ridgeplot_data_prep <- function(data, axis_offset = NULL) {
    value <- NULL # NSE vars
    data <- data[, -1] # remove cell_ID col
    data <- data.table::melt(data,
        measure.vars = colnames(data),
        variable.name = "feat",
        value.name = "value"
    )
    # apply axis offset
    if (!is.null(axis_offset)) {
        data[, value := value + axis_offset]
    }
    data
}

# return a `gg` with ridgeplot
# data must be a long format table with two columns ("feat" and "value")
.ridgeplot_single <- function(data,
    scale_axis = "identity",
    ridge_alpha = 1,
    ridge_scale = 5,
    stat = "density_ridges",
    remove_na = NULL,
    background_color = "white",
    title = "",
    show_legend = FALSE,
    theme_param = list()
) {
    # accumulate plot aesthetics
    plot_aes <- aes_string2(
        x = "value",
        y = "feat",
        fill = "feat"
    )

    # accumulate theme param
    theme_param$panel.background <- theme_param$panel.background %null%
        ggplot2::element_rect(fill = background_color)

    ggplot2::ggplot(data = data) +
        ggridges::geom_density_ridges(
            mapping = plot_aes,
            alpha = ridge_alpha,
            stat = stat,
            na.rm = remove_na %null% FALSE,
            scale = ridge_scale,
            show.legend = show_legend
        ) +
        ggplot2::labs(x = "Value", y = "Feature", title = title) +
        ggplot2::scale_x_continuous(transform = scale_axis) +
        ggridges::theme_ridges() +
        do.call(ggplot2::theme, theme_param)
}


