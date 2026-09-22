#' @title showClusterHeatmap
#' @name showClusterHeatmap
#' @description Creates heatmap based on identified clusters
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values expression values to use
#' (e.g. "normalized", "scaled", "custom")
#' @param feats vector of features to use, default to 'all'
#' @param cluster_column name of column to use for clusters
#' (e.g. "leiden_clus")
#' @param cor correlation score to calculate distance
#' (e.g. "pearson", "spearman")
#' @param distance distance method to use for hierarchical clustering,
#' default to "ward.D"
#' @inheritParams gmulti_params
#' @inheritDotParams ComplexHeatmap::Heatmap
#' @details Correlation heatmap of selected clusters.
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' showClusterHeatmap(g, cluster_column = "leiden_clus")
#'
#' @export
showClusterHeatmap <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    feats = "all",
    cluster_column,
    cor = c("pearson", "spearman"),
    distance = "ward.D",
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "showClusterHeatmap",
    view = NULL,
    space = NULL,
    ...) {
    # package Check
    package_check(pkg_name = "ComplexHeatmap", repository = "Bioc")
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


    ## correlation
    cor <- match.arg(cor, c("pearson", "spearman"))
    values <- match.arg(expression_values, c("normalized", "scaled", "custom"))

    ## subset expression data
    if (feats[1] != "all") {
        all_feats <- gobject@feat_ID[[feat_type]]
        detected_feats <- feats[feats %in% all_feats]
    } else {
        # NULL = all feats in calculateMetaTable()
        detected_feats <- NULL
    }

    metatable <- calculateMetaTable(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        metadata_cols = cluster_column,
        selected_feats = detected_feats
    )
    dcast_metatable <- data.table::dcast.data.table(metatable,
        formula = variable ~ uniq_ID,
        value.var = "value"
    )
    testmatrix <- dt_to_matrix(x = dcast_metatable)

    # correlation
    cormatrix <- cor_flex(x = testmatrix, method = cor)
    cordist <- stats::as.dist(1 - cormatrix, diag = TRUE, upper = TRUE)
    corclus <- stats::hclust(d = cordist, method = distance)

    hmap <- ComplexHeatmap::Heatmap(
        matrix = cormatrix,
        cluster_rows = corclus,
        cluster_columns = corclus,
        ...
    )

    return(plot_output_handler(
        gobject = gobject,
        plot_object = hmap,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}








#' @title plotHeatmap
#' @name plotHeatmap
#' @description Creates heatmap for features and clusters.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_params
#' @param expression_values expression values to use
#' (e.g. "normalized", "scaled", "custom")
#' @param feats features to use
#' @param cluster_column name of column to use for clusters
#' (e.g. "leiden_clus")
#' @param cluster_order method to determine cluster order
#' (e.g. "size", "correlation", "custom")
#' @param cluster_custom_order custom order for clusters
#' @param cluster_color_code color code for clusters
#' @param cluster_cor_method method for cluster correlation,
#' default to "pearson"
#' @param cluster_hclust_method method for hierarchical clustering of
#' clusters, default to "ward.D"
#' @param feat_order method to determine features order
#'  (e.g. "correlation", "custom")
#' @param feat_custom_order custom order for features
#' @param feat_cor_method method for features correlation, default to "pearson"
#' @param feat_hclust_method method for hierarchical clustering of
#' features, default to "complete"
#' @param show_values which values to show on heatmap
#' (e.g. "rescaled", "z-scaled", "original")
#' @param size_vertical_lines sizes for vertical lines
#' @param gradient_colors deprecated
#' @param feat_label_selection subset of features to show on y-axis
#' @param axis_text_y_size size for y-axis text
#' @param legend_nrows number of rows for the cluster legend
#' @inheritParams gmulti_params
#' @returns ggplot
#' @details If you want to display many features there are 2 ways to proceed:
#' \itemize{
#'   \item{1. set axis_text_y_size to a very small value and show all features}
#'   \item{2. provide a subset of features to display to feat_label_selection}
#' }
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotHeatmap(g,
#'     feats = c("Gm19935", "Gna12", "Ccnd2", "Btbd17"),
#'     cluster_column = "leiden_clus"
#' )
#'
#' @export
plotHeatmap <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    feats,
    cluster_column = NULL,
    cluster_order = c("size", "correlation", "custom"),
    cluster_custom_order = NULL,
    cluster_color_code = NULL,
    cluster_cor_method = "pearson",
    cluster_hclust_method = "ward.D",
    feat_order = c("correlation", "custom"),
    feat_custom_order = NULL,
    feat_cor_method = "pearson",
    feat_hclust_method = "complete",
    show_values = c("rescaled", "z-scaled", "original"),
    size_vertical_lines = 1.1,
    gradient_colors = deprecated(),
    gradient_color = NULL,
    gradient_style = c("divergent", "sequential"),
    feat_label_selection = NULL,
    axis_text_y_size = NULL,
    legend_nrows = 1,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "plotHeatmap",
    view = NULL,
    space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression", "spatial_enrichment"))
    # deprecate
    if (GiottoUtils::is_present(gradient_colors)) {
        deprecate_warn(
            "0.0.0.9000",
            "GiottoVisuals::plotHeatmap(gradient_colors = )",
            "GiottoVisuals::plotHeatmap(gradient_color = )"
        )
        gradient_color <- gradient_colors
    }

    show_values <- match.arg(show_values,
        choices = c("rescaled", "z-scaled", "original")
    )

    heatmap_data <- .create_heatmap_dt(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = expression_values,
        feats = feats,
        cluster_column = cluster_column,
        cluster_order = cluster_order,
        cluster_custom_order = cluster_custom_order,
        cluster_cor_method = cluster_cor_method,
        cluster_hclust_method = cluster_hclust_method,
        feat_order = feat_order,
        feat_custom_order = feat_custom_order,
        feat_cor_method = feat_cor_method,
        feat_hclust_method = feat_hclust_method
    )

    cell_order_DT <- heatmap_data[["cell_DT"]]
    subset_values_DT <- heatmap_data[["DT"]]
    x_lines <- heatmap_data[["x_lines"]]

    ## assign colors to each cluster
    if (is.null(cluster_color_code)) {
        clus_values <- unique(cell_order_DT[[cluster_column]])
        clus_colors <- set_default_color_discrete_heatmap_clus(
            instrs = instructions(gobject))(n = length(clus_values))
        names(clus_colors) <- clus_values
    } else {
        clus_colors <-
            set_default_color_discrete(
                colors = cluster_color_code,
                NULL,
                NULL,
                "interpolate"
            )(n = length(unique(cell_order_DT[[cluster_column]])))
    }

    ## bar on top ##
    clus_pl <- ggplot2::ggplot()
    clus_pl <- clus_pl +
        ggplot2::geom_raster(
            data = cell_order_DT,
            aes_string2(
                x = "cells",
                y = 1,
                fill = cluster_column))
    clus_pl <- clus_pl +
        ggplot2::geom_vline(
            xintercept = x_lines,
            color = "white",
            linewidth = size_vertical_lines
        )
    clus_pl <- clus_pl +
        ggplot2::scale_fill_manual(
            values = clus_colors,
            guide = ggplot2::guide_legend(title = "", nrow = legend_nrows)
        )
    clus_pl <- clus_pl + ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        legend.position = "top",
        plot.margin = ggplot2::margin(0, 0, 0, 5.5, "pt")
    )
    clus_pl <- clus_pl + ggplot2::labs(x = "", y = "clusters")

    ## rescale values ##
    if (show_values == "original") {
        value_column <- "expression"
        midpoint <- max(subset_values_DT$expression) / 2
    } else if (show_values == "z-scaled") {
        value_column <- "z_scores"
        midpoint <- max(subset_values_DT$z_scores) / 2
    } else {
        value_column <- "scale_scores"
        midpoint <- 0.5
    }

    # create empty plot
    empty <- ggplot2::ggplot()
    empty <- empty + ggplot2::theme_classic()
    empty <- empty + ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 0, 0, 0, "cm")
    )

    ### heatmap ###
    hmap <- ggplot2::ggplot()
    hmap <- hmap +
        ggplot2::geom_tile(
            data = subset_values_DT,
            aes_string2(
                x = "cells",
                y = "feats",
                fill = value_column
            )
        )
    hmap <- hmap +
        ggplot2::geom_vline(
            xintercept = x_lines,
            color = "white",
            linewidth = size_vertical_lines
        )
    hmap <- hmap + set_default_color_continuous_heatmap(
        colors = gradient_color,
        instrs = instructions(gobject),
        midpoint = midpoint,
        style = gradient_style,
        guide = ggplot2::guide_colorbar(title = ""),
        type = "fill"
    )

    if (is.null(feat_label_selection)) {
        hmap <- hmap + ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = axis_text_y_size),
            plot.margin = ggplot2::margin(0, 0, 5.5, 5.5, "pt")
        )

        ## align and combine
        aligned <- cowplot::align_plots(
            clus_pl, empty, hmap +
            ggplot2::theme(legend.position = "none"),
            align = "v", axis = "l"
        )
        aligned <- append(aligned, list(cowplot::get_legend(hmap)))
        combplot <- cowplot::plot_grid(
            plotlist = aligned,
            ncol = 2, rel_widths = c(1, 0.2),
            nrow = 2, rel_heights = c(0.2, 1)
        )
    } else {
        # set defaults
        if (is.null(axis_text_y_size)) {
            axis_text_y_size <- 0.8 * 11 / ggplot2::.pt
        }

        # finish heatmap
        hmap <- hmap + ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 5.5, 5.5, "pt")
        )
        ### axis ###
        featDT <- subset_values_DT[, c("feats"), with = FALSE]
        featDT <- unique(data.table::setorder(featDT, feats))

        # data.table variables
        featOrder <- subset_feats <- NULL

        featDT[, featOrder := seq_len(.N)]
        featDT[, subset_feats := ifelse(feats %in% feat_label_selection,
            as.character(feats), ""
        )]

        axis <- ggplot2::ggplot(
            data = featDT,
            aes(x = 0, y = featOrder,label = subset_feats))
        axis <- axis + ggrepel::geom_text_repel(
            min.segment.length = grid::unit(0, "pt"),
            color = "grey30", ## ggplot2 theme_grey() axis text
            size = axis_text_y_size ## default ggplot2 theme_grey() axis text
        )
        axis <- axis + ggplot2::scale_x_continuous(
            limits = c(0, 1), expand = c(0, 0),
            breaks = NULL, labels = NULL, name = NULL
        )
        axis <- axis + ggplot2::scale_y_continuous(
            limits = c(0, nrow(featDT)), expand = c(0, 0),
            breaks = NULL, labels = NULL, name = NULL
        )
        axis <- axis + ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(0, 0, 0, 0, "cm")
        )

        ## align and combine
        aligned <- cowplot::align_plots(
            clus_pl, empty, empty, hmap +
            theme(legend.position = "none"),
            axis,
            align = "h",
            axis = "b"
        )
        aligned <- append(aligned, list(cowplot::get_legend(hmap)))
        combplot <- cowplot::plot_grid(
            plotlist = aligned,
            ncol = 3, rel_widths = c(1, 0.2, 0.1),
            nrow = 2, rel_heights = c(0.2, 1)
        )
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = combplot,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}






#' @title plotMetaDataHeatmap
#' @name plotMetaDataHeatmap
#' @description Creates heatmap for features within aggregated clusters.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_params
#' @param expression_values expression values to use
#' (e.g. "normalized", "scaled", "custom")
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param selected_feats subset of features to use
#' @param first_meta_col if more than 1 metadata column,
#' select the x-axis factor
#' @param second_meta_col if more than 1 metadata column,
#' select the facetting factor
#' @param show_values which values to show on heatmap
#' (e.g. "zscores", "original", "zscores_rescaled")
#' @param custom_cluster_order custom cluster order (default = NULL)
#' @param clus_cor_method correlation method for clusters, default to "pearson"
#' @param clus_cluster_method hierarchical cluster method for the clusters,
#' default to "complete"
#' @param custom_feat_order custom feature order (default = NULL)
#' @param feat_cor_method correlation method for features, default to "pearson"
#' @param feat_cluster_method hierarchical cluster method for the features,
#' default to "complete"
#' @param x_text_size size of x-axis text
#' @param x_text_angle angle of x-axis text
#' @param y_text_size size of y-axis text
#' @param strip_text_size size of strip text
#' @param plot_title deprecated. Use title param
#' @inheritParams gmulti_params
#' @details Creates heatmap for the average expression of selected features
#' in the different annotation/cluster groups.
#' Calculation of cluster or feature order is done on the provided expression
#' values, but visualization
#' is by default on the z-scores. Other options are the original values or
#' z-scores rescaled per feature (-1 to 1).
#' @seealso \code{\link{plotMetaDataCellsHeatmap}} for numeric cell annotation
#' instead of feature expression.
#' @returns ggplot or data.table
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotMetaDataHeatmap(g,
#'     metadata_cols = "leiden_clus",
#'     selected_feats = c("Gna12", "Ccnd2", "Btbd17", "Gm19935")
#' )
#'
#' @export
plotMetaDataHeatmap <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        metadata_cols = NULL,
        selected_feats = NULL,
        first_meta_col = NULL,
        second_meta_col = NULL,
        show_values = c("zscores", "original", "zscores_rescaled"),
        custom_cluster_order = NULL,
        clus_cor_method = "pearson",
        clus_cluster_method = "complete",
        custom_feat_order = NULL,
        feat_cor_method = "pearson",
        feat_cluster_method = "complete",
        gradient_color = NULL,
        gradient_midpoint = 0,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        x_text_size = 10,
        x_text_angle = 45,
        y_text_size = 10,
        strip_text_size = 8,
        title = NULL,
        plot_title = deprecated(),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "plotMetaDataHeatmap",
        view = NULL,
        space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression", "spatial_enrichment"))
    # deprecate
    if (GiottoUtils::is_present(plot_title)) {
        deprecate_warn(
            "0.0.0.9000",
            "GiottoVisuals::plotMetaDataHeatmap(plot_title = )",
            "GiottoVisuals::plotMetaDataHeatmap(title = )"
        )
        title <- plot_title
    }

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

    metaDT <- calculateMetaTable(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = expression_values,
        metadata_cols = metadata_cols,
        selected_feats = selected_feats
    )

    # data.table variables
    zscores <- value <- NULL

    metaDT[, "zscores" := scale(value), by = c("variable")]
    metaDT[, "zscores_rescaled_per_feat" := scales::rescale(zscores,
        to = c(-1, 1)
    ), by = c("variable")]

    show_values <- match.arg(show_values,
        choices = c("zscores", "original", "zscores_rescaled")
    )
    if (show_values == "zscores") {
        show_values <- "zscores"
    } else if (show_values == "original") {
        show_values <- "value"
    } else {
        show_values <- "zscores_rescaled_per_feat"
    }

    ## visualization
    if (length(metadata_cols) > 2) {
        message("visualization is only possible for 1 or 2 metadata annotations,
        data.table is returned")
        return(metaDT)
    }


    ## order of feats and clusters

    main_factor <- ifelse(length(metadata_cols) == 1L,
        metadata_cols, first_meta_col
    )
    testmain <- metaDT[, mean(value), by = c("variable", main_factor)]
    testmain_matrix <- dt_dcast_string(testmain, "variable", main_factor, "V1")
    testmain_mat <- as.matrix(testmain_matrix[, -1])
    rownames(testmain_mat) <- testmain_matrix$variable

    # for clusters
    if (is.null(custom_cluster_order)) {
        cormatrix <- cor_flex(x = testmain_mat, method = clus_cor_method)
        cordist <- stats::as.dist(1 - cormatrix, diag = TRUE, upper = TRUE)
        corclus <- stats::hclust(d = cordist, method = clus_cluster_method)
        clus_names <- rownames(cormatrix)
        names(clus_names) <- seq_len(length(clus_names))
        clus_sort_names <- clus_names[corclus$order]
    } else {
        clus_sort_names <- unique(as.character(custom_cluster_order))
        if (all(colnames(testmain_mat) %in% clus_sort_names) == FALSE) {
            stop("\n custom cluster order is given,
                but not all clusters are represented \n")
        }
    }


    # for feats
    if (is.null(custom_feat_order)) {
        feat_cormatrix <- cor_flex(
            x = t(testmain_mat),
            method = feat_cor_method
        )
        feat_cordist <- stats::as.dist(1 - feat_cormatrix,
            diag = TRUE, upper = TRUE
        )
        feat_corclus <- stats::hclust(
            d = feat_cordist,
            method = feat_cluster_method
        )
        feat_names <- rownames(feat_cormatrix)
        names(feat_names) <- seq_len(length(feat_names))
        feat_sort_names <- feat_names[feat_corclus$order]
    } else {
        feat_sort_names <- unique(as.character(custom_feat_order))
        if (all(rownames(testmain_mat) %in% feat_sort_names) == FALSE) {
            stop("\n custom feat order is given,
                but not all feats are represented \n")
        }
    }

    if (length(metadata_cols) == 1) {
        # data.table variables
        factor_column <- variable <- NULL

        metaDT[, factor_column := factor(get(metadata_cols),
            levels = clus_sort_names
        )]
        metaDT[, variable := factor(get("variable"),
            levels = feat_sort_names
        )]

        # set gradient information
        if (!is.null(gradient_limits) &
            is.vector(gradient_limits) & length(gradient_limits) == 2) {
            lower_lim <- gradient_limits[[1]]
            upper_lim <- gradient_limits[[2]]

            numeric_data <- metaDT[[show_values]]
            limit_numeric_data <- ifelse(numeric_data > upper_lim, upper_lim,
                ifelse(numeric_data < lower_lim, lower_lim, numeric_data)
            )
            metaDT[[show_values]] <- limit_numeric_data
        }


        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::geom_tile(
            data = metaDT,
            aes_string2(
                x = "factor_column",
                y = "variable",
                fill = show_values
            ),
            color = "black"
        )
        pl <- pl + set_default_color_continuous_heatmap(
            colors = gradient_color,
            instrs = instructions(gobject),
            midpoint = gradient_midpoint,
            style = gradient_style,
            type = "fill"
        )
        pl <- pl + ggplot2::theme_classic()
        pl <- pl + ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                size = x_text_size,
                angle = x_text_angle,
                hjust = 1, vjust = 1
            ),
            axis.text.y = ggplot2::element_text(size = y_text_size),
            legend.title = ggplot2::element_blank()
        )
        pl <- pl + ggplot2::labs(x = metadata_cols, y = "feats")

        if (!is.null(title)) {
            pl <- pl + ggplot2::ggtitle(title) +
                ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        }


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
                c(list(
                    gobject = gobject, plot_object = pl,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(pl)
        }
    } else {
        if (is.null(first_meta_col) | is.null(second_meta_col)) {
            message("both first_meta_col and second_meta_col need to be defined,
            return data.table")
            return(metaDT)
        } else {
            # data.table variables
            factor_1_column <- factor_2_column <- variable <- NULL

            metaDT[, factor_1_column := factor(
                get(first_meta_col),
                clus_sort_names
            )]
            metaDT[, factor_2_column := as.factor(get(second_meta_col))]
            metaDT[, variable := factor(get("variable"),
                levels = feat_sort_names
            )]

            # set gradient information
            if (!is.null(gradient_limits) &
                is.vector(gradient_limits) & length(gradient_limits) == 2) {
                lower_lim <- gradient_limits[[1]]
                upper_lim <- gradient_limits[[2]]

                numeric_data <- metaDT[[show_values]]
                limit_numeric_data <- ifelse(numeric_data > upper_lim,
                    upper_lim,
                    ifelse(numeric_data < lower_lim, lower_lim, numeric_data)
                )
                metaDT[[show_values]] <- limit_numeric_data
            }

            pl <- ggplot2::ggplot()
            pl <- pl +
                ggplot2::geom_tile(
                    data = metaDT,
                    aes_string2(
                        x = "factor_1_column",
                        y = "variable",
                        fill = show_values
                    ),
                    color = "black"
                )
            pl <- pl + set_default_color_continuous_heatmap(
                colors = gradient_color,
                instrs = instructions(gobject),
                midpoint = gradient_midpoint,
                style = gradient_style,
                type = "fill"
            )
            pl <- pl +
                ggplot2::facet_grid(stats::reformulate("factor_2_column"))
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    size = x_text_size,
                    angle = x_text_angle,
                    hjust = 1, vjust = 1
                ),
                axis.text.y = ggplot2::element_text(size = y_text_size),
                strip.text = ggplot2::element_text(size = strip_text_size),
                legend.title = ggplot2::element_blank()
            )
            pl <- pl +
                ggplot2::labs(
                    x = first_meta_col,
                    y = "feats",
                    title = second_meta_col
                )


            # print, return and save parameters
            show_plot <- ifelse(is.na(show_plot),
                readGiottoInstructions(
                    gobject,
                    param = "show_plot"
                ),
                show_plot
            )
            save_plot <- ifelse(is.na(save_plot),
                readGiottoInstructions(
                    gobject,
                    param = "save_plot"
                ),
                save_plot
            )
            return_plot <- ifelse(is.na(return_plot),
                readGiottoInstructions(
                    gobject,
                    param = "return_plot"
                ),
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
                    c(list(
                        gobject = gobject,
                        plot_object = pl,
                        default_save_name = default_save_name
                    ), save_param)
                )
            }

            ## return plot
            if (return_plot == TRUE) {
                return(pl)
            }
        }
    }
}


#' @title plotMetaDataCellsHeatmap
#' @name plotMetaDataCellsHeatmap
#' @description Creates heatmap for numeric cell metadata within
#' aggregated clusters.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_params
#' @param midpoint deprecated. Use gradient_midpoint.
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param spat_enr_names spatial enrichment results to include
#' @param value_cols value columns to use
#' @param first_meta_col if more than 1 metadata column,
#' select the x-axis factor
#' @param second_meta_col if more than 1 metadata column,
#' select the facetting factor
#' @param show_values which values to show on heatmap
#' (e.g. "zscores", "original", "zscores_rescaled")
#' @param custom_cluster_order custom cluster order (default = NULL)
#' @param clus_cor_method correlation method for clusters, default to "pearson"
#' @param clus_cluster_method hierarchical cluster method for the clusters,
#' default to "complete"
#' @param custom_values_order custom values order (default = NULL)
#' @param values_cor_method correlation method for values, default to "pearson"
#' @param values_cluster_method hierarchical cluster method for the values,
#' default to "complete"
#' @param x_text_size size of x-axis text
#' @param x_text_angle angle of x-axis text
#' @param y_text_size size of y-axis text
#' @param y_axis_title character. default is "genes".
#' @param strip_text_size size of strip text
#' @inheritParams gmulti_params
#' @details Creates heatmap for the average values of selected value columns
#' in the different annotation groups.
#' @seealso \code{\link{plotMetaDataHeatmap}} for feature expression instead
#' of numeric cell annotation data.
#' @returns ggplot or data.table
#'
#' @export
plotMetaDataCellsHeatmap <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    metadata_cols = NULL,
    spat_enr_names = NULL,
    value_cols = NULL,
    first_meta_col = NULL,
    second_meta_col = NULL,
    show_values = c("zscores", "original", "zscores_rescaled"),
    custom_cluster_order = NULL,
    clus_cor_method = "pearson",
    clus_cluster_method = "complete",
    custom_values_order = NULL,
    values_cor_method = "pearson",
    values_cluster_method = "complete",
    gradient_color = NULL,
    gradient_midpoint = 0,
    gradient_style = c("divergent", "sequential"),
    midpoint = deprecated(),
    x_text_size = 8,
    x_text_angle = 45,
    y_text_size = 8,
    y_axis_title = "genes",
    strip_text_size = 8,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "plotMetaDataCellsHeatmap",
    view = NULL,
    space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression", "spatial_enrichment"))
    # deprecate
    if (GiottoUtils::is_present(midpoint)) {
        deprecate_warn(
            "0.0.0.9000",
            "GiottoVisuals::plotMetaDataCellsHeatmap(midpoint = )",
            "GiottoVisuals::plotMetaDataCellsHeatmap(gradient_midpoint = )"
        )
        gradient_midpoint <- midpoint
    }

    metaDT <- calculateMetaTableCells(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        value_cols = value_cols,
        metadata_cols = metadata_cols,
        spat_enr_names = spat_enr_names
    )

    # data.table variables
    zscores <- value <- NULL

    metaDT[, "zscores" := scale(value), by = c("variable")]
    metaDT[, "zscores_rescaled_per_gene" := scales::rescale(zscores,
        to = c(-1, 1)
    ),
    by = c("variable")
    ]

    show_values <- match.arg(show_values,
        choices = c("zscores", "original", "zscores_rescaled")
    )
    if (show_values == "zscores") {
        show_values <- "zscores"
    } else if (show_values == "original") {
        show_values <- "value"
    } else {
        show_values <- "zscores_rescaled_per_gene"
    }

    ## visualization
    if (length(metadata_cols) > 2) {
        message("visualization is only possible for 1 or 2 metadata annotations,
        data.table is returned")
        return(metaDT)
    }


    ## order of genes and clusters

    main_factor <- ifelse(length(metadata_cols) == 1,
        metadata_cols, first_meta_col
    )
    testmain <- metaDT[, mean(value), by = c("variable", main_factor)]
    testmain_matrix <- dt_dcast_string(testmain, "variable", main_factor, "V1")
    testmain_mat <- as.matrix(testmain_matrix[, -1])
    rownames(testmain_mat) <- testmain_matrix$variable

    # for clusters
    if (is.null(custom_cluster_order)) {
        cormatrix <- cor_flex(x = testmain_mat, method = clus_cor_method)
        cordist <- stats::as.dist(1 - cormatrix,
            diag = TRUE, upper = TRUE
        )
        corclus <- stats::hclust(d = cordist, method = clus_cluster_method)
        clus_names <- rownames(cormatrix)
        names(clus_names) <- seq_len(length(clus_names))
        clus_sort_names <- clus_names[corclus$order]
    } else {
        clus_sort_names <- unique(as.character(custom_cluster_order))
        if (all(colnames(testmain_mat) %in% clus_sort_names) == FALSE) {
            stop("\n custom cluster order is given, but not all clusters
                are represented \n")
        }
    }


    # for genes
    if (is.null(custom_values_order)) {
        values_cormatrix <- cor_flex(
            x = t(testmain_mat),
            method = values_cor_method
        )
        values_cordist <- stats::as.dist(1 - values_cormatrix,
            diag = TRUE,
            upper = TRUE
        )
        values_corclus <- stats::hclust(
            d = values_cordist,
            method = values_cluster_method
        )
        values_names <- rownames(values_cormatrix)
        names(values_names) <- seq_len(length(values_names))
        values_sort_names <- values_names[values_corclus$order]
    } else {
        values_sort_names <- unique(as.character(custom_values_order))
        if (all(rownames(testmain_mat) %in% values_sort_names) == FALSE) {
            stop("\n custom values order is given,
                but not all values are represented \n")
        }
    }

    if (length(metadata_cols) == 1) {
        # data.table variables
        factor_column <- variable <- NULL

        metaDT[, factor_column := factor(get(metadata_cols),
            levels = clus_sort_names
        )]
        metaDT[, variable := factor(get("variable"),
            levels = values_sort_names
        )]

        pl <- ggplot2::ggplot()
        pl <- pl +
            ggplot2::geom_tile(
                data = metaDT,
                aes_string2(
                    x = "factor_column",
                    y = "variable",
                    fill = show_values
                ),
                color = "black"
            )
        pl <- pl + set_default_color_continuous_heatmap(
            colors = gradient_color,
            instrs = instructions(gobject),
            midpoint = gradient_midpoint,
            style = gradient_style,
            type = "fill"
        )
        pl <- pl + ggplot2::theme_classic()
        pl <- pl + ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                size = x_text_size,
                angle = x_text_angle,
                hjust = 1, vjust = 1
            ),
            axis.text.y = ggplot2::element_text(size = y_text_size),
            legend.title = ggplot2::element_blank()
        )
        pl <- pl + ggplot2::labs(x = metadata_cols, y = y_axis_title)
    } else {
        if (is.null(first_meta_col) | is.null(second_meta_col)) {
            message("both first_meta_col and second_meta_col need to be defined,
            return data.table")
            return(metaDT)
        } else {
            # data.table variables
            factor_1_column <- factor_2_column <- variable <- NULL

            metaDT[, factor_1_column := factor(
                get(first_meta_col),
                clus_sort_names
            )]
            metaDT[, factor_2_column := as.factor(get(second_meta_col))]
            metaDT[, variable := factor(get("variable"),
                levels = values_sort_names
            )]

            pl <- ggplot2::ggplot()
            pl <- pl +
                ggplot2::geom_tile(
                    data = metaDT,
                    aes_string2(
                        x = "factor_1_column",
                        y = "variable",
                        fill = show_values),
                    color = "black"
                )
            pl <- pl + set_default_color_continuous_heatmap(
                colors = gradient_color,
                instrs = instructions(gobject),
                midpoint = gradient_midpoint,
                style = gradient_style,
                type = "fill"
            )
            pl <- pl +
                ggplot2::facet_grid(stats::reformulate("factor_2_column"))
            pl <- pl + ggplot2::theme_classic()
            pl <- pl + ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                    size = x_text_size,
                    angle = x_text_angle,
                    hjust = 1, vjust = 1
                ),
                axis.text.y = ggplot2::element_text(size = y_text_size),
                strip.text = ggplot2::element_text(size = strip_text_size),
                legend.title = ggplot2::element_blank()
            )
            pl <- pl + ggplot2::labs(
                x = first_meta_col, y = y_axis_title,
                title = second_meta_col
            )
        }
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




# helpers ####



#' @title Create heatmap data.table
#' @name .create_heatmap_dt
#' @description creates order for clusters
#' @inheritParams data_access_params
#' @param expression_values expression values to use
#' (e.g. "normalized", "scaled", "custom")
#' @param feats features to use
#' @param cluster_column name of column to use for clusters
#' (e.g. "leiden_clus")
#' @param cluster_order method to determine cluster order
#' (e.g. "size", "correlation", "custom")
#' @param cluster_custom_order custom order for clusters
#' @param cluster_cor_method method for cluster correlation,
#' default to "pearson"
#' @param cluster_hclust_method method for hierarchical clustering of clusters,
#' default to "ward.D"
#' @param feat_order method to determine features order
#' (e.g. "correlation", "custom")
#' @param feat_custom_order custom order for features
#' @param feat_cor_method method for features correlation, default to "pearson"
#' @param feat_hclust_method method for hierarchical clustering of features,
#' default to "complete"
#' @return list
#' @details Creates input data.tables for plotHeatmap function.
#' @keywords internal
.create_heatmap_dt <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        feats,
        cluster_column = NULL,
        cluster_order = c("size", "correlation", "custom"),
        cluster_custom_order = NULL,
        cluster_cor_method = "pearson",
        cluster_hclust_method = "ward.D",
        feat_order = c("correlation", "custom"),
        feat_custom_order = NULL,
        feat_cor_method = "pearson",
        feat_hclust_method = "complete") {
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


    # epxression data
    values <- match.arg(
        expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )
    expr_values <- getExpression(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        values = values,
        output = "matrix"
    )

    # subset expression data
    detected_feats <- feats[feats %in% rownames(expr_values)]
    subset_values <- expr_values[rownames(expr_values) %in% detected_feats, ]

    # metadata
    cell_metadata <- pDataDT(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    # feat order
    feat_order <- match.arg(feat_order, c("correlation", "custom"))


    ### cluster order ###
    clus_sort_names <- .decide_cluster_order(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = expression_values,
        feats = feats,
        cluster_column = cluster_column,
        cluster_order = cluster_order,
        cluster_custom_order = cluster_custom_order,
        cor_method = cluster_cor_method,
        hclust_method = cluster_hclust_method
    )

    ## data.table ##
    # Coerce through dgCMatrix rather than calling `as.matrix()` on whatever
    # `getExpression()` returned. A disk-backed store has no `as.matrix()`
    # method, and giving it one would be the wrong fix: `as.matrix()` means
    # dense, so a backed carrier would densify with none of the size checks
    # its own materialization path applies. Going via dgCMatrix keeps that
    # path, and every in-memory carrier used here -- base matrix, dgCMatrix,
    # DelayedArray, IterableMatrix -- answers it too.
    #
    # `melt_matrix()` needs a base array (it calls `as.data.frame.table()`),
    # so the result is dense either way; only the route differs.
    subset_values_DT <- melt_matrix(
        as.matrix(methods::as(subset_values, "dgCMatrix")),
        varnames = c("feats", "cells"),
        value.name = "expression",
        as.is = TRUE
    )
    subset_values_DT <- merge(
        subset_values_DT,
        by.x = "cells",
        cell_metadata[, c("cell_ID", cluster_column), with = FALSE],
        by.y = "cell_ID"
    )
    subset_values_DT[, (cluster_column) := factor(
        x = get(cluster_column),
        levels = clus_sort_names
    )]

    # NSE vars
    z_scores <- scale_scores <- V1 <- cells <- NULL

    subset_values_DT[, feats := factor(feats, unique(detected_feats))]
    subset_values_DT[, z_scores := scale(expression), by = feats]
    subset_values_DT[, scale_scores := scales::rescale(
        x = expression,
        to = c(0, 1)
    ),
    by = feats
    ]


    ## order cells by mean expression ##
    cell_order_DT <- subset_values_DT[, mean(expression),
        by = c("cells", cluster_column)
    ]
    cell_order_DT <- cell_order_DT[order(get(cluster_column), V1)]
    subset_values_DT[, cells := factor(cells, cell_order_DT$cells)]

    ## get x-coordines for vertical lines in heatmap
    x_lines <- cumsum(as.vector(table(cell_order_DT[[cluster_column]])))

    ## order feats ##
    if (feat_order == "correlation") {
        featsum_per_clus <- subset_values_DT[, sum(expression),
            by = c("feats", cluster_column)
        ]

        my_formula <- paste0("feats~", cluster_column)
        test_mat <- data.table::dcast.data.table(
            data = featsum_per_clus,
            formula = my_formula,
            value.var = "V1"
        )
        test_matrix <- as.matrix(test_mat[, -1])
        rownames(test_matrix) <- test_mat$feats

        feat_dist <- stats::as.dist(1 - cor_flex(t_flex(test_matrix),
            method = feat_cor_method
        ))
        feat_clus <- stats::hclust(feat_dist, method = feat_hclust_method)

        feat_labels <- rownames(test_matrix)
        feat_index <- seq_len(length(feat_labels))
        names(feat_index) <- feat_labels

        final_feat_order <- names(feat_index[match(
            feat_clus$order,
            feat_index
        )])
        subset_values_DT[, "feats" := factor(feats, final_feat_order)]
    } else if (feat_order == "custom") {
        if (is.null(feat_custom_order)) {
            .gstop("with custom feat order,
                the feat_custom_order parameter needs to be provided")
        }
        subset_values_DT[, "feats" := factor(feats, feat_custom_order)]
    }

    cell_order_DT[["cells"]] <- factor(cell_order_DT[["cells"]],
        levels = as.character(cell_order_DT[["cells"]])
    )

    return(
        list(
            DT = subset_values_DT,
            x_lines = x_lines,
            cell_DT = cell_order_DT
        )
    )
}
