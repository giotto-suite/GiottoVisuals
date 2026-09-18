#' @title violinPlot
#' @name violinPlot
#' @description Creates violinplot for selected clusters
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values expression values to
#' use (e.g. "normalized", "scaled", "custom")
#' @param feats features to plot
#' @param cluster_column name of column to use for
#' clusters (e.g. "leiden_clus")
#' @param color_violin color violin according to "genes" or "clusters"
#' @param cluster_custom_order custom order of clusters
#' @param cluster_color_code color code for clusters
#' @param strip_position position of gene
#' labels (e.g. "top", "right", "left", "bottom")
#' @param strip_text size of strip text
#' @param axis_text_x_size size of x-axis text
#' @param axis_text_y_size size of y-axis text
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' violinPlot(g,
#'     feats = c("Gna12", "Ccnd2", "Btbd17"),
#'     cluster_column = "leiden_clus"
#' )
#'
#' @export
violinPlot <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    feats = NULL,
    cluster_column,
    cluster_custom_order = NULL,
    color_violin = c("feats", "cluster"),
    cluster_color_code = NULL,
    strip_position = c("top", "right", "left", "bottom"),
    strip_text = 7,
    axis_text_x_size = 10,
    axis_text_y_size = 6,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "violinPlot",
    view = NULL,
    space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression", "spatial_enrichment"))
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

    ## strip position
    strip_position <- match.arg(
        strip_position,
        c("top", "right", "left", "bottom")
    )

    ## color of violin plots
    color_violin <- match.arg(color_violin, c("feats", "cluster"))

    ## expression data ##
    values <- match.arg(
        expression_values,
        unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
    )

    # collect expression values using spatValues
    expr_dt <- spatValues(gobject,
        feats = unique(feats),
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        verbose = FALSE
    )
    selected_feats <- setdiff(colnames(expr_dt), "cell_ID")

    # stop if no feats have been found
    if (length(selected_feats) == 0) {
        stop("No overlapping features have been found,
            check inputer for parameter 'feats'")
    }

    # cluster column
    clus_dt <- spatValues(gobject,
        feats = cluster_column,
        spat_unit = spat_unit,
        feat_type = feat_type,
        verbose = FALSE
    )

    metadata_expr <- clus_dt[expr_dt, on = "cell_ID"]


    metadata_expr_m <-
        data.table::melt.data.table(metadata_expr,
            measure.vars = selected_feats,
            variable.name = "feats"
        )
    metadata_expr_m[, feats := factor(feats, selected_feats)]
    metadata_expr_m[[cluster_column]] <-
        as.factor(metadata_expr_m[[cluster_column]])

    if (!is.null(cluster_custom_order)) {
        cluster_custom_order <- unique(cluster_custom_order)
        metadata_expr_m[[cluster_column]] <-
            factor(
                x = metadata_expr_m[[cluster_column]],
                levels = cluster_custom_order
            )
    }


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    if (color_violin == "feats") {
        pl <- pl + ggplot2::geom_violin(
            data = metadata_expr_m,
            aes_string2(
                x = cluster_column,
                y = "value", 
                fill = "feats"
            ),
            width = 1, 
            scale = "width",
            show.legend = FALSE
        )
    } else {
        pl <- pl + ggplot2::geom_violin(
            data = metadata_expr_m,
            aes_string2(
                x = cluster_column, 
                y = "value",
                fill = cluster_column
            ),
            width = 1, 
            scale = "width",
            show.legend = FALSE
        )

        # provide own color scheme for clusters
        if (!is.null(cluster_color_code)) {
            pl <- pl + ggplot2::scale_fill_manual(values = cluster_color_code)
        }
    }


    pl <- pl + ggplot2::facet_wrap(. ~ feats,
        ncol = 1,
        strip.position = strip_position
    )
    pl <- pl + ggplot2::theme(
        strip.text = element_text(size = strip_text),
        axis.text.x = element_text(
            size = axis_text_x_size,
            angle = 45, hjust = 1, vjust = 1
        ),
        axis.text.y = element_text(size = axis_text_y_size)
    )
    pl <- pl + ggplot2::labs(x = "", y = "normalized expression")



    # print, return and save parameters
    show_plot <- ifelse(is.null(show_plot),
        readGiottoInstructions(gobject, param = "show_plot"),
        show_plot
    )
    save_plot <- ifelse(is.null(save_plot),
        readGiottoInstructions(gobject, param = "save_plot"),
        save_plot
    )
    return_plot <- ifelse(is.null(return_plot),
        readGiottoInstructions(gobject, param = "return_plot"),
        return_plot
    )

    ## print plot
    if (show_plot == TRUE) {
        print(pl)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call(
            "all_plots_save_function",
            c(
                list(
                    gobject = gobject,
                    plot_object = pl, default_save_name = default_save_name
                ),
                save_param
            )
        )
    }

    ## return plot
    if (return_plot == TRUE) {
        return(pl)
    }
}
