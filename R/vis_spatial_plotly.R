# * ####
## 3-D plotly ####
## ----------- ##

# ** dimension plot ####


#' @title .dimPlot_2d_plotly
#' @name .dimPlot_2d_plotly
#' @description Visualize cells at their 2D dimension reduction coordinates
#' with plotly
#' @returns plotly object
#' @keywords internal
.dimPlot_2d_plotly <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        spat_enr_names = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 0.5,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        color_as_factor = TRUE,
        cell_color = NULL,
        cell_color_code = NULL,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        label_size = 4,
        edge_alpha = NULL,
        point_size = 5) {
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

    # data.table variables
    cell_ID <- NULL

    ## dimension reduction ##
    dim_dfr <- getDimReduction(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = "cells",
        reduction_method = dim_reduction_to_use,
        name = dim_reduction_name,
        output = "data.table"
    )

    dim_dfr <- dim_dfr[, c(dim1_to_use, dim2_to_use)]
    dim_names <- colnames(dim_dfr)
    dim_DT <- data.table::as.data.table(dim_dfr)
    dim_DT[, "cell_ID" := rownames(dim_dfr)]


    ## annotated cell metadata
    cell_metadata <- combineMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spat_enr_names = spat_enr_names
    )
    annotated_DT <- merge(cell_metadata, dim_DT, by = "cell_ID")


    # create input for network
    if (show_NN_network == TRUE) {
        # nn_network
        selected_nn_network <- getNearestNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            nn_type = nn_network_to_use,
            name = network_name,
            output = "igraph"
        )
        network_DT <- data.table::as.data.table(igraph::as_data_frame(
            selected_nn_network,
            what = "edges"
        ))

        # annotated network
        old_dim_names <- dim_names

        annotated_network_DT <- merge(network_DT, dim_DT,
            by.x = "from",
            by.y = "cell_ID"
        )
        from_dim_names <- paste0("from_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = from_dim_names
        )

        annotated_network_DT <- merge(annotated_network_DT, dim_DT,
            by.x = "to", by.y = "cell_ID"
        )
        to_dim_names <- paste0("to_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = to_dim_names
        )
    }


    if (dim_reduction_to_use == "pca") {
        pca_object <- getDimReduction(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "dimObj"
        )
        eigenvalues <- slot(pca_object, "misc")$eigenvalues

        if (!is.null(eigenvalues)) {
            total <- sum(eigenvalues)
            var_expl_vec <- (eigenvalues / total) * 100
            dim1_x_variance <- var_expl_vec[dim1_to_use]
            dim2_y_variance <- var_expl_vec[dim2_to_use]
        }
    }


    if (!is.null(select_cells) & !is.null(select_cell_groups)) {
        if (is.null(cell_color)) {
            stop("\n selection of cells is based on cell_color paramter,
                which is a metadata column \n")
        }
        message("You have selected both individual cell IDs and a group
        of cells")
        group_cell_IDs <- annotated_DT[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- annotated_DT[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
    }


    if (!is.null(select_cells)) {
        annotated_DT_other <- annotated_DT[!annotated_DT$cell_ID %in%
            select_cells]
        annotated_DT_selected <- annotated_DT[annotated_DT$cell_ID %in%
            select_cells]

        if (show_NN_network == TRUE) {
            annotated_network_DT <- annotated_network_DT[
                annotated_network_DT$to %in% select_cells &
                    annotated_network_DT$from %in% select_cells
            ]
        }

        # if specific cells are selected
        # annotated_DT = annotated_DT_selected
    }


    ## if no subsets are required
    if (is.null(select_cells) & is.null(select_cell_groups)) {
        annotated_DT_selected <- annotated_DT
        annotated_DT_other <- NULL
    }


    ## annotated_DT_selected = all selected cells or all cells if no selection
    ## annotated_DT_other = all not selected cells or NULL if no selection


    pl <- plotly::plot_ly()
    if (show_NN_network == TRUE) {
        if (is.null(edge_alpha)) {
            edge_alpha <- 0.5
        } else if (is.character(edge_alpha)) {
            warning("Edge_alpha for plotly mode is not adjustable yet.
                    Default 0.5 will be set\n")
            edge_alpha <- 0.5
        }

        pl <- pl %>% plotly::add_segments(
            name = network_name,
            type = "scatter",
            x = annotated_network_DT[[from_dim_names[1]]],
            y = annotated_network_DT[[from_dim_names[2]]],
            xend = annotated_network_DT[[to_dim_names[1]]],
            yend = annotated_network_DT[[to_dim_names[2]]],
            line = list(
                color = "lightgray",
                width = 0.5
            ),
            opacity = edge_alpha
        )
    }

    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl %>% plotly::add_trace(
            type = "scatter", mode = "markers",
            x = annotated_DT_selected[[dim_names[1]]],
            y = annotated_DT_selected[[dim_names[2]]],
            color = cell_color,
            colors = cell_color,
            marker = list(size = point_size)
        )
    } else if (cell_color %in% colnames(annotated_DT_selected)) {
        if (is.null(cell_color_code)) {
            number_colors <- length(unique(annotated_DT[[cell_color]]))
            cell_color_code <- set_default_color_discrete_cell(
                instrs = instructions(gobject)
            )(n = number_colors)
        }
        if (color_as_factor) {
            annotated_DT_selected[[cell_color]] <- as.factor(
                annotated_DT_selected[[cell_color]]
            )
        }


        pl <- pl %>% plotly::add_trace(
            type = "scatter", mode = "markers",
            x = annotated_DT_selected[[dim_names[1]]],
            y = annotated_DT_selected[[dim_names[2]]],
            color = annotated_DT_selected[[cell_color]],
            colors = cell_color_code,
            legendgroup = annotated_DT_selected[[cell_color]],
            marker = list(size = point_size)
        )

        if (!is.null(select_cells) & show_other_cells) {
            pl <- pl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT_other[[dim_names[1]]],
                y = annotated_DT_other[[dim_names[2]]],
                # legendgroup = annotated_DT[[cell_color]],
                marker = list(
                    size = other_point_size,
                    color = other_cell_color
                ),
                showlegend = FALSE
            )
        }

        if (show_cluster_center == TRUE | show_center_label == TRUE) {
            annotated_DT_centers <- annotated_DT_selected[, .(
                center_1 = stats::median(get(dim_names[1])),
                center_2 = stats::median(get(dim_names[2]))
            ),
            by = cell_color
            ]
            annotated_DT_centers[[cell_color]] <- as.factor(
                annotated_DT_centers[[cell_color]]
            )
            if (show_cluster_center == TRUE) {
                pl <- pl %>% plotly::add_trace(
                    type = "scatter", mode = "markers",
                    x = annotated_DT_centers[["center_1"]],
                    y = annotated_DT_centers[["center_2"]],
                    color = annotated_DT_centers[[cell_color]],
                    colors = cell_color_code,
                    legendgroup = annotated_DT_centers[[cell_color]],
                    marker = list(
                        size = center_point_size, symbol = "x",
                        symbols = "x"
                    ),
                    showlegend = FALSE
                )
            }

            if (show_center_label == TRUE) {
                pl <- pl %>%
                    plotly::add_text(
                        x = annotated_DT_centers[["center_1"]],
                        y = annotated_DT_centers[["center_2"]],
                        type = "scatter", mode = "text",
                        text = annotated_DT_centers[[cell_color]],
                        textposition = "middle right",
                        textfont = list(color = "#000000", size = 16),
                        showlegend = FALSE
                    )
            }
        }
    } else {
        stop("cell_color does not exist!\n")
    }



    if (dim_reduction_to_use == "pca") {
        if (!is.null(eigenvalues)) {
            x_name <- paste0("pca", "-", dim_names[1])
            y_name <- paste0("pca", "-", dim_names[2])
            x_title <- sprintf(
                "%s explains %.02f%% of variance",
                x_name, var_expl_vec[1]
            )
            y_title <- sprintf(
                "%s explains %.02f%% of variance", y_name,
                var_expl_vec[2]
            )
        }
    } else {
        x_title <- paste(dim_reduction_to_use, dim_names[1], sep = " ")
        y_title <- paste(dim_reduction_to_use, dim_names[2], sep = " ")
    }
    pl <- pl %>% plotly::layout(
        xaxis = list(title = x_title),
        yaxis = list(title = y_title),
        legend = list(x = 100, y = 0.5, font = list(
            family = "sans-serif",
            size = 12
        ))
    )

    return(pl)
}


#' @title .dimPlot_3d_plotly
#' @name .dimPlot_3d_plotly
#' @description Visualize cells at their 3D dimension reduction coordinates
#' with plotly
#' @returns plotly object
#' @keywords internal
.dimPlot_3d_plotly <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        dim3_to_use = 3,
        spat_enr_names = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 0.5,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        color_as_factor = TRUE,
        cell_color = NULL,
        cell_color_code = NULL,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        label_size = 4,
        edge_alpha = NULL,
        point_size = 1) {
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

    # data.table variables
    cell_ID <- NULL

    ## dimension reduction ##
    dim_mat <- getDimReduction(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = "cells",
        reduction_method = dim_reduction_to_use,
        name = dim_reduction_name,
        output = "matrix"
    )
    dim_mat <- dim_mat[, c(dim1_to_use, dim2_to_use, dim3_to_use)]
    dim_names <- colnames(dim_mat)
    dim_DT <- data.table::as.data.table(dim_mat, keep.rownames = TRUE)
    data.table::setnames(dim_DT, old = "rn", new = "cell_ID")


    ## annotated cell metadata
    cell_metadata <- combineMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spat_enr_names = spat_enr_names
    )
    annotated_DT <- merge(cell_metadata, dim_DT, by = "cell_ID")


    # create input for network
    if (show_NN_network == TRUE) {
        # nn_network
        selected_nn_network <- getNearestNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            nn_type = nn_network_to_use,
            name = network_name,
            output = "igraph"
        )
        network_DT <- data.table::as.data.table(igraph::as_data_frame(
            selected_nn_network,
            what = "edges"
        ))

        # annotated network
        old_dim_names <- dim_names

        annotated_network_DT <- merge(network_DT, dim_DT,
            by.x = "from",
            by.y = "cell_ID"
        )
        from_dim_names <- paste0("from_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = from_dim_names
        )

        annotated_network_DT <- merge(annotated_network_DT, dim_DT,
            by.x = "to", by.y = "cell_ID"
        )
        to_dim_names <- paste0("to_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = to_dim_names
        )
    }

    if (dim_reduction_to_use == "pca") {
        pca_object <- getDimReduction(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "dimObj"
        )

        eigenvalues <- slot(pca_object, "misc")$eigenvalues
        if (!is.null(eigenvalues)) {
            total <- sum(eigenvalues)
            var_expl_vec <- (eigenvalues / total) * 100
            dim1_x_variance <- var_expl_vec[dim1_to_use]
            dim2_y_variance <- var_expl_vec[dim2_to_use]
        }
    }

    ## create subsets if needed
    if (!is.null(select_cells) & !is.null(select_cell_groups)) {
        if (is.null(cell_color)) {
            stop("\n selection of cells is based on cell_color parameter,
                which is a metadata column \n")
        }
        message("You have selected both individual cell IDs and a group of
        cells")
        group_cell_IDs <- annotated_DT[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- annotated_DT[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
    }

    if (!is.null(select_cells)) {
        annotated_DT_other <- annotated_DT[!annotated_DT$cell_ID %in%
            select_cells]
        annotated_DT_selected <- annotated_DT[annotated_DT$cell_ID %in%
            select_cells]

        if (show_NN_network == TRUE) {
            annotated_network_DT <- annotated_network_DT[
                annotated_network_DT$to %in% select_cells &
                    annotated_network_DT$from %in% select_cells
            ]
        }

        # if specific cells are selected
        annotated_DT <- annotated_DT_selected
    }

    ## if no subsets are required
    if (is.null(select_cells) & is.null(select_cell_groups)) {
        annotated_DT_selected <- annotated_DT
        annotated_DT_other <- NULL
    }

    ## annotated_DT_selected = all selected cells or all cells if no selection
    ## annotated_DT_other = all not selected cells or NULL if no selection


    pl <- plotly::plot_ly()
    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl %>% plotly::add_trace(
            type = "scatter3d", mode = "markers",
            x = annotated_DT_selected[[dim_names[1]]],
            y = annotated_DT_selected[[dim_names[2]]],
            z = annotated_DT_selected[[dim_names[3]]],
            color = cell_color,
            colors = cell_color,
            marker = list(size = 2),
            legendgroup = annotated_DT_selected[[cell_color]]
        )
    } else {
        if (cell_color %in% colnames(annotated_DT_selected)) {
            if (is.null(cell_color_code)) {
                number_colors <- length(
                    unique(annotated_DT_selected[[cell_color]])
                )
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instructions(gobject)
                )(n = number_colors)
            }
            if (color_as_factor) {
                annotated_DT_selected[[cell_color]] <- as.factor(
                    annotated_DT_selected[[cell_color]]
                )
            }

            pl <- pl %>% plotly::add_trace(
                type = "scatter3d", mode = "markers",
                x = annotated_DT_selected[[dim_names[1]]],
                y = annotated_DT_selected[[dim_names[2]]],
                z = annotated_DT_selected[[dim_names[3]]],
                color = annotated_DT_selected[[cell_color]],
                colors = cell_color_code,
                marker = list(size = point_size),
                legendgroup = annotated_DT_selected[[cell_color]]
            )

            if (!is.null(select_cells) & show_other_cells) {
                pl <- pl %>% plotly::add_trace(
                    type = "scatter3d", mode = "markers",
                    x = annotated_DT_other[[dim_names[1]]],
                    y = annotated_DT_other[[dim_names[2]]],
                    z = annotated_DT_other[[dim_names[3]]],
                    # colors = other_cell_color,
                    marker = list(
                        size = other_point_size,
                        color = other_cell_color
                    ),
                    showlegend = FALSE
                )
            }


            if (show_cluster_center == TRUE | show_center_label == TRUE) {
                annotated_DT_centers <- annotated_DT_selected[, .(
                    center_1 = stats::median(get(dim_names[1])),
                    center_2 = stats::median(get(dim_names[2])),
                    center_3 = stats::median(get(dim_names[3]))
                ),
                by = cell_color
                ]
                annotated_DT_centers[[cell_color]] <- as.factor(
                    annotated_DT_centers[[cell_color]]
                )
                if (show_cluster_center == TRUE) {
                    pl <- pl %>% plotly::add_trace(
                        mode = "markers",
                        type = "scatter3d",
                        data = annotated_DT_centers,
                        x = ~center_1,
                        y = ~center_2,
                        z = ~center_3,
                        color = annotated_DT_centers[[cell_color]],
                        colors = cell_color_code,
                        inherit = FALSE,
                        marker = list(size = 2, symbol = "x", symbols = "x"),
                        legendgroup = annotated_DT_centers[[cell_color]],
                        showlegend = FALSE
                    )
                }
                if (show_center_label == TRUE) {
                    message(" center label is not clear to see in 3D plot\n You
                    can shut it down with show_center_label = FALSE")
                    pl <- pl %>% plotly::add_trace(
                        mode = "text",
                        type = "scatter3d",
                        data = annotated_DT_centers,
                        x = ~center_1,
                        y = ~center_2,
                        z = ~center_3,
                        text = annotated_DT_centers[[cell_color]],
                        legendgroup = annotated_DT_centers[[cell_color]],
                        inherit = FALSE,
                        showlegend = FALSE
                    )
                }
            }
        } else {
            stop("cell_color does not exist!\n")
        }
    }

    if (show_NN_network) {
        edges <- plotly_network(
            annotated_network_DT,
            "from_Dim.1", "from_Dim.2", "from_Dim.3",
            "to_Dim.1", "to_Dim.2", "to_Dim.3"
        )
        if (is.null(edge_alpha)) {
            edge_alpha <- 0.5
        } else if (is.character(edge_alpha)) {
            warning("Edge_alpha for plotly mode is not adjustable yet.
                    Default 0.5 will be set\n")
            edge_alpha <- 0.5
        }

        pl <- pl %>% plotly::add_trace(
            name = network_name,
            mode = "lines",
            type = "scatter3d",
            data = edges,
            x = ~x, y = ~y, z = ~z,
            inherit = FALSE,
            line = list(color = "lightgray", width = 0.5),
            opacity = edge_alpha
        )
    }

    if (dim_reduction_to_use == "pca") {
        if (!is.null(eigenvalues)) {
            x_name <- paste0("pca", "-", dim_names[1])
            y_name <- paste0("pca", "-", dim_names[2])
            z_name <- paste0("pca", "-", dim_names[3])
            x_title <- sprintf(
                "%s explains %.02f%% of variance",
                x_name, var_expl_vec[1]
            )
            y_title <- sprintf(
                "%s explains %.02f%% of variance",
                y_name, var_expl_vec[2]
            )
            z_title <- sprintf(
                "%s explains %.02f%% of variance",
                z_name, var_expl_vec[3]
            )
        }
    } else {
        x_title <- paste(dim_reduction_to_use, dim_names[1], sep = " ")
        y_title <- paste(dim_reduction_to_use, dim_names[2], sep = " ")
        z_title <- paste(dim_reduction_to_use, dim_names[3], sep = " ")
    }
    pl <- pl %>% plotly::layout(
        scene = list(
            xaxis = list(title = x_title),
            yaxis = list(title = y_title),
            zaxis = list(title = z_title)
        ),
        legend = list(
            x = 100, y = 0.5,
            font = list(family = "sans-serif", size = 12)
        )
    )
    return(pl)
}







#' @rdname dimPlot
#' @returns plotly (dimplot3D only)
#' @export
dimPlot3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        dim3_to_use = 3,
        spat_enr_names = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 2,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        color_as_factor = TRUE,
        cell_color = NULL,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        label_size = 4,
        edge_alpha = NULL,
        point_size = 3,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "dim3D") {
    if (is.null(dim3_to_use)) {
        message("create 2D plot")

        pl <- .dimPlot_2d_plotly(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            dim_reduction_to_use = dim_reduction_to_use,
            dim_reduction_name = dim_reduction_name,
            dim1_to_use = dim1_to_use,
            dim2_to_use = dim2_to_use,
            spat_enr_names = spat_enr_names,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_NN_network = show_NN_network,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            color_as_factor = color_as_factor,
            cell_color = cell_color,
            cell_color_code = cell_color_code,
            show_cluster_center = show_cluster_center,
            show_center_label = show_center_label,
            center_point_size = center_point_size,
            label_size = label_size,
            edge_alpha = edge_alpha,
            point_size = point_size
        )
    } else {
        message("create 3D plot")
        pl <- .dimPlot_3d_plotly(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            dim_reduction_to_use = dim_reduction_to_use,
            dim_reduction_name = dim_reduction_name,
            dim1_to_use = dim1_to_use,
            dim2_to_use = dim2_to_use,
            dim3_to_use = dim3_to_use,
            spat_enr_names = spat_enr_names,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_NN_network = show_NN_network,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            color_as_factor = color_as_factor,
            cell_color = cell_color,
            cell_color_code = cell_color_code,
            show_cluster_center = show_cluster_center,
            show_center_label = show_center_label,
            center_point_size = center_point_size,
            label_size = label_size,
            edge_alpha = edge_alpha,
            point_size = point_size
        )
    }


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
}


#' @title plotUMAP_3D
#' @name plotUMAP_3D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_name name of UMAP
#' @param default_save_name default save name of UMAP plot
#' @inheritDotParams dimPlot3D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters.
#' @family reduced dimension visualizations
#' @returns plotly
#' @examples
#' g <- GiottoData::loadGiottoMini("starmap")
#' plotUMAP_3D(g, dim_reduction_name = "3D_umap")
#'
#' @export
plotUMAP_3D <- function(
        gobject,
        dim_reduction_name = "umap",
        default_save_name = "UMAP_3D",
        ...) {
    dimPlot3D(
        gobject = gobject,
        dim_reduction_to_use = "umap",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}


#' @title plotTSNE_3D
#' @name plotTSNE_3D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_name name of TSNE
#' @param default_save_name default save name of TSNE plot
#' @inheritDotParams dimPlot3D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters.
#' @family reduced dimension visualizations
#' @returns plotly
#' @examples
#' g <- GiottoData::loadGiottoMini("starmap")
#' plotTSNE_3D(g)
#'
#' @export
plotTSNE_3D <- function(
        gobject,
        dim_reduction_name = "tsne",
        default_save_name = "TSNE_3D",
        ...) {
    dimPlot3D(
        gobject = gobject,
        dim_reduction_to_use = "tsne",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}


#' @title plotPCA_3D
#' @name plotPCA_3D
#' @description Visualize cells according to 3D PCA dimension reduction
#' @param gobject giotto object
#' @param dim_reduction_name name of PCA
#' @param default_save_name default save name of PCA plot
#' @inheritDotParams dimPlot3D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters.
#' @family reduced dimension visualizations
#' @returns plotly
#' @examples
#' g <- GiottoData::loadGiottoMini("starmap")
#' plotPCA_3D(g)
#'
#' @export
plotPCA_3D <- function(
        gobject,
        dim_reduction_name = "pca",
        default_save_name = "PCA_3D",
        ...) {
    dimPlot3D(
        gobject = gobject,
        dim_reduction_to_use = "pca",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}






# ** ####
# ** spatial 3D plot ####

#' @title .spatPlot_2d_plotly
#' @name .spatPlot_2d_plotly
#' @description Visualize cells at their 2D spatial locations with plotly
#' @returns plotly object
#' @keywords internal
.spatPlot_2d_plotly <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spat_loc_name = "raw",
        sdimx = NULL,
        sdimy = NULL,
        spat_enr_names = NULL,
        point_size = 3,
        cell_color = NULL,
        cell_color_code = NULL,
        color_as_factor = TRUE,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 0.5,
        show_network = FALSE,
        spatial_network_name = "spatial_network",
        network_color = "lightgray",
        network_alpha = 1,
        other_cell_alpha = 0.5,
        show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        grid_color = NULL,
        grid_alpha = 1,
        show_legend = TRUE,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        show_plot = FALSE) {
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

    ## get spatial cell locations
    cell_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table"
    )
    if (is.null(cell_locations)) {
        return(NULL)
    }


    ## extract spatial network
    if (show_network) {
        # Edges plus their endpoint coordinates. The network holds no
        # geometry itself, so the positions are attached here from the
        # locations -- and because that join is inner, a narrowed object
        # drops the edges whose cells are gone with no work here.
        spatial_network <- GiottoClass::annotateSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            spatial_network_name = spatial_network_name,
            spat_loc_name = spat_loc_name,
            coordinates = TRUE
        )
    } else {
        spatial_network <- NULL
    }

    ## extract spatial grid
    if (show_grid == TRUE) {
        spatial_grid <- getSpatialGrid(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = spatial_grid_name
        )
    } else {
        spatial_grid <- NULL
    }

    ## get cell metadata
    cell_metadata <- try(
        expr = combineMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            spat_loc_name = spat_loc_name,
            spat_enr_names = spat_enr_names
        ),
        silent = TRUE
    )


    if (inherits(cell_metadata, "try-error")) {
        cell_locations_metadata <- cell_locations
    } else if (nrow(cell_metadata) == 0) {
        cell_locations_metadata <- cell_locations
    } else {
        cell_locations_metadata <- cell_metadata
    }

    ## create subsets if needed
    if (!is.null(select_cells) & !is.null(select_cell_groups)) {
        message("You have selected both individual cell IDs and a group
        of cells")
        group_cell_IDs <- cell_locations_metadata[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- cell_locations_metadata[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
    }


    if (!is.null(select_cells)) {
        cell_locations_metadata_other <- cell_locations_metadata[
            !cell_locations_metadata$cell_ID %in% select_cells
        ]
        cell_locations_metadata_selected <- cell_locations_metadata[
            cell_locations_metadata$cell_ID %in% select_cells
        ]
        spatial_network <- spatial_network[spatial_network$to %in%
            select_cells & spatial_network$from %in%
            select_cells]

        # if specific cells are selected
        # cell_locations_metadata = cell_locations_metadata_selected
    } else if (is.null(select_cells)) {
        cell_locations_metadata_selected <- cell_locations_metadata
        cell_locations_metadata_other <- NULL
    }



    ### set scale
    axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))

    ### set ratio
    ratio <- plotly_axis_scale_2D(cell_locations,
        sdimx = sdimx,
        sdimy = sdimy,
        mode = axis_scale,
        custom_ratio = custom_ratio
    )



    pl <- plotly::plot_ly()

    ## create network
    if (show_network == TRUE) {
        if (is.null(spatial_network)) {
            stop("No usable spatial network specified! Please choose a
                network with spatial_network_name=xxx")
        } else {
            if (is.null(network_alpha)) {
                network_alpha <- 0.5
            } else if (is.character(network_alpha)) {
                warning("Edge_alpha for plotly mode is not adjustable yet.
                        Default 0.5 will be set\n")
                network_alpha <- 0.5
            }
            pl <- pl %>% plotly::add_segments(
                name = spatial_network_name,
                type = "scatter",
                x = spatial_network[["sdimx_begin"]],
                y = spatial_network[["sdimy_begin"]],
                xend = spatial_network[["sdimx_end"]],
                yend = spatial_network[["sdimy_end"]],
                line = list(
                    color = network_color,
                    width = 0.5
                ),
                opacity = network_alpha
            )
        }
    }

    ## create grid
    if (show_grid == TRUE) {
        if (is.null(spatial_grid)) {
            stop("No usable spatial grid specified! Please choose a
                network with spatial_grid_name=xxx")
        } else {
            if (is.null(grid_color)) {
                grid_color <- "black"
            }
            edges <- plotly_grid(spatial_grid)
            pl <- pl %>% plotly::add_segments(
                name = "spatial_grid",
                type = "scatter",
                data = edges,
                x = ~x,
                y = ~y,
                xend = ~x_end,
                yend = ~y_end,
                line = list(
                    color = grid_color,
                    width = 1
                ),
                opacity = grid_alpha
            )
        }
    }



    if (!is.null(cell_color)) {
        if (cell_color %in% colnames(cell_locations_metadata_selected)) {
            if (is.null(cell_color_code)) {
                number_colors <- length(unique(
                    cell_locations_metadata_selected[[cell_color]]
                ))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instructions(gobject)
                )(n = number_colors)
            }
            cell_locations_metadata_selected[[cell_color]] <- as.factor(
                cell_locations_metadata_selected[[cell_color]]
            )
            pl <- pl %>% plotly::add_trace(
                type = "scatter",
                mode = "markers",
                x = cell_locations_metadata_selected[[sdimx]],
                y = cell_locations_metadata_selected[[sdimy]],
                color = cell_locations_metadata_selected[[cell_color]],
                colors = cell_color_code,
                marker = list(size = point_size)
            )


            if (!is.null(select_cells) & show_other_cells) {
                pl <- pl %>% plotly::add_trace(
                    type = "scatter",
                    mode = "markers",
                    data = cell_locations_metadata_other,
                    name = "unselected cells",
                    x = ~sdimx,
                    y = ~sdimy,
                    marker = list(
                        size = other_point_size,
                        color = other_cell_color
                    ),
                    opacity = other_cell_alpha
                )
            }
        } else {
            message("cell_color does not exist!")
        }
    } else {
        pl <- pl %>% plotly::add_trace(
            type = "scatter",
            mode = "markers",
            name = "selected cells",
            x = cell_locations_metadata_selected[[sdimx]],
            y = cell_locations_metadata_selected[[sdimy]],
            colors = "lightblue",
            marker = list(size = point_size)
        )

        if (!is.null(select_cells) & show_other_cells) {
            pl <- pl %>% plotly::add_trace(
                type = "scatter",
                mode = "markers",
                data = cell_locations_metadata_other,
                name = "unselected cells",
                x = ~sdimx,
                y = ~sdimy,
                marker = list(
                    size = other_point_size,
                    color = other_cell_color
                ),
                opacity = other_cell_alpha
            )
        }
    }


    pl <- pl %>%
        plotly::layout(
            list(
                xaxis = list(title = "X", nticks = x_ticks),
                yaxis = list(title = "Y", nticks = y_ticks)
            ),
            legend = list(
                x = 100, y = 0.5,
                font = list(
                    family = "sans-serif",
                    size = 12
                )
            )
        )


    return((pl))
}



#' @title .spatPlot_3d_plotly
#' @name .spatPlot_3d_plotly
#' @description Visualize cells at their 3D spatial locations with plotly
#' @returns plotly object
#' @keywords internal
.spatPlot_3d_plotly <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spat_loc_name = "raw",
        sdimx = NULL,
        sdimy = NULL,
        sdimz = NULL,
        spat_enr_names = NULL,
        point_size = 3,
        point_alpha = 1,
        cell_color = NULL,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        color_as_factor = TRUE,
        gradient_limits = NULL,
        gradient_style = "divergent",
        gradient_midpoint = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 3,
        show_network = FALSE,
        spatial_network_name = "spatial_network",
        network_color = NULL,
        network_alpha = 1,
        other_cell_alpha = 0.5,
        show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        title = "",
        show_legend = TRUE,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        show_plot = FALSE) {
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

    ## get spatial cell locations
    cell_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table"
    )
    if (is.null(cell_locations)) {
        return(NULL)
    }

    ## extract spatial network
    if (show_network) {
        # Edges plus their endpoint coordinates. The network holds no
        # geometry itself, so the positions are attached here from the
        # locations -- and because that join is inner, a narrowed object
        # drops the edges whose cells are gone with no work here.
        spatial_network <- GiottoClass::annotateSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            spatial_network_name = spatial_network_name,
            spat_loc_name = spat_loc_name,
            coordinates = TRUE
        )
    } else {
        spatial_network <- NULL
    }

    ## extract spatial grid
    if (show_grid) {
        spatial_grid <- getSpatialGrid(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = spatial_grid_name
        )
    } else {
        spatial_grid <- NULL
    }

    ## get cell metadata
    cell_metadata <- try(
        expr = combineMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            spat_loc_name = spat_loc_name,
            spat_enr_names = spat_enr_names
        ),
        silent = TRUE
    )


    if (inherits(cell_metadata, "try-error")) {
        cell_locations_metadata <- cell_locations
    } else if (nrow(cell_metadata) == 0) {
        cell_locations_metadata <- cell_locations
    } else {
        cell_locations_metadata <- cell_metadata
    }


    ## create subsets if needed
    if (!is.null(select_cells) && !is.null(select_cell_groups)) {
        message("You have selected both individual cell IDs and a group of
        cells")
        group_cell_IDs <- cell_locations_metadata[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- cell_locations_metadata[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
    }

    if (!is.null(select_cells)) {
        cell_locations_metadata_other <- cell_locations_metadata[
            !cell_locations_metadata$cell_ID %in% select_cells
        ]
        cell_locations_metadata_selected <- cell_locations_metadata[
            cell_locations_metadata$cell_ID %in% select_cells
        ]
        spatial_network <- spatial_network[spatial_network$to %in%
            select_cells & spatial_network$from %in% select_cells]

        # if specific cells are selected
        # cell_locations_metadata = cell_locations_metadata_selected
    } else if (is.null(select_cells)) {
        cell_locations_metadata_selected <- cell_locations_metadata
        cell_locations_metadata_other <- NULL
    }



    ### set scale
    axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))

    ### set ratio
    ratio <- plotly_axis_scale_3D(cell_locations,
        sdimx = sdimx,
        sdimy = sdimy,
        sdimz = sdimz,
        mode = axis_scale,
        custom_ratio = custom_ratio
    )



    pl <- plotly::plot_ly()

    ## plot points
    pl <- giotto_point_3d(pl,
        data = cell_locations_metadata_selected,
        cell_color = cell_color,
        color_as_factor = color_as_factor,
        cell_color_code = cell_color_code,
        cell_color_gradient = cell_color_gradient,
        gradient_limits = gradient_limits,
        gradient_style = gradient_style,
        gradient_midpoint = gradient_midpoint,
        point_size = point_size,
        point_alpha = point_alpha,
        data_other = cell_locations_metadata_other,
        select_cells = select_cells,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        other_point_size = other_point_size,
        other_cell_alpha = other_cell_alpha,
        instrs = instructions(gobject)
    )


    ## plot spatial network
    if (!is.null(spatial_network) & show_network == TRUE) {
        if (is.null(network_color)) {
            network_color <- "red"
        }
        edges <- plotly_network(spatial_network)

        pl <- pl %>% plotly::add_trace(
            name = "sptial network",
            mode = "lines",
            type = "scatter3d",
            data = edges,
            x = ~x,
            y = ~y,
            z = ~z,
            line = list(color = network_color, width = 0.5),
            opacity = network_alpha
        )
    }

    ## plot spatial grid
    # 3D grid is not clear to view


    pl <- pl %>%
        plotly::layout(
            scene = list(
                xaxis = list(title = "X", nticks = x_ticks),
                yaxis = list(title = "Y", nticks = y_ticks),
                zaxis = list(title = "Z", nticks = z_ticks),
                aspectmode = "manual",
                aspectratio = list(
                    x = ratio[[1]],
                    y = ratio[[2]],
                    z = ratio[[3]]
                )
            ),
            legend = list(
                x = 100, y = 0.5,
                font = list(family = "sans-serif", size = 12)
            )
        )


    return(pl)
}





#' @rdname spatPlot
#' @param sdimz z-axis dimension name (default = 'sdimy'). Can be set to `NULL`
#' to force 2D plotting with {plotly}
#' @param grid_alpha opacity of spatial grid
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#' @export
spatPlot3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        sdimz = "sdimz",
        spat_enr_names = NULL,
        point_size = 3,
        point_alpha = 1,
        cell_color = NULL,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        color_as_factor = TRUE,
        gradient_limits = NULL,
        gradient_style = "divergent",
        gradient_midpoint = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 3,
        other_cell_alpha = 0.5,
        show_network = FALSE,
        spatial_network_name = "Delaunay_network",
        network_color = NULL,
        network_alpha = 1,
        show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        grid_color = NULL,
        grid_alpha = 1,
        title = "",
        show_legend = TRUE,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spat3D") {
    if (is.null(sdimz)) {
        vmsg(.is_debug = TRUE, "create 2D plot")

        pl <- .spatPlot_2d_plotly(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sdimx = sdimx,
            sdimy = sdimy,
            spat_enr_names = spat_enr_names,
            point_size = point_size,
            cell_color = cell_color,
            color_as_factor = color_as_factor,
            cell_color_code = cell_color_code,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_network = show_network,
            network_color = network_color,
            network_alpha = network_alpha,
            other_cell_alpha = other_cell_alpha,
            spatial_network_name = spatial_network_name,
            show_grid = show_grid,
            grid_color = grid_color,
            grid_alpha = grid_alpha,
            spatial_grid_name = spatial_grid_name,
            show_legend = show_legend,
            axis_scale = axis_scale,
            custom_ratio = custom_ratio,
            x_ticks = x_ticks,
            y_ticks = y_ticks,
            show_plot = FALSE
        )
    } else {
        vmsg(.is_debug = TRUE, "create 3D plot")
        pl <- .spatPlot_3d_plotly(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            sdimx = sdimx,
            sdimy = sdimy,
            sdimz = sdimz,
            spat_enr_names = spat_enr_names,
            point_size = point_size,
            point_alpha = point_alpha,
            cell_color = cell_color,
            color_as_factor = color_as_factor,
            cell_color_code = cell_color_code,
            cell_color_gradient = cell_color_gradient,
            gradient_limits = gradient_limits,
            gradient_style = gradient_style,
            gradient_midpoint = gradient_midpoint,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_network = show_network,
            network_color = network_color,
            network_alpha = network_alpha,
            other_cell_alpha = other_cell_alpha,
            spatial_network_name = spatial_network_name,
            spatial_grid_name = spatial_grid_name,
            show_legend = show_legend,
            axis_scale = axis_scale,
            custom_ratio = custom_ratio,
            x_ticks = x_ticks,
            y_ticks = y_ticks,
            z_ticks = z_ticks,
            show_plot = FALSE
        )
    }

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
}









# ** ####
# ** spatial & dimension 3D plot ####

#' @title spatDimPlot3D
#' @name spatDimPlot3D
#' @description Visualize cells according to spatial AND dimension
#' reduction coordinates in plotly mode
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#'
#' @param spat_loc_name name for spatial locations
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param sdimz = spatial dimension to use on z-axis
#'
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param show_cluster_center show the center of each cluster
#' @param show_center_label provide a label for each cluster
#' @param center_point_size size of the center point
#' @param label_size size of the center label
#'
#' @param select_cell_groups select subset of cells/clusters based on
#' cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param dim_point_size size of points in dim. reduction space
#' @param nn_network_color color of nn network
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#'
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_grid_alpha alpha of spatial grid
#' @param spatial_point_size size of spatial points
#' @param spatial_network_color color of spatial network
#' @param spatial_network_alpha alpha of spatial network
#'
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#' @param legend_text_size size of legend
#' @returns plotly
#' @details Description of parameters.
#' @family spatial and dimension reduction visualizations
#' @export
spatDimPlot3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        plot_alignment = c("horizontal", "vertical"),
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        dim3_to_use = 3,
        spat_loc_name = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        sdimz = "sdimz",
        spat_enr_names = NULL,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        nn_network_color = "lightgray",
        nn_network_alpha = 0.5,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        label_size = 16,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1.5,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        dim_point_size = 3,
        show_spatial_network = FALSE,
        spatial_network_name = "Delaunay_network",
        spatial_network_color = "lightgray",
        spatial_network_alpha = 0.5,
        show_spatial_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        spatial_grid_color = NULL,
        spatial_grid_alpha = 0.5,
        spatial_point_size = 3,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        legend_text_size = 12,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatDimPlot3D") {
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

    # data.table variables
    cell_ID <- NULL

    plot_alignment <- match.arg(plot_alignment,
        choices = c("horizontal", "vertical")
    )

    # ********data prepare********#
    ## dimension reduction ##
    dim_dfr <- getDimReduction(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = "cells",
        reduction_method = dim_reduction_to_use,
        name = dim_reduction_name,
        output = "matrix"
    )
    dim_dfr <- dim_dfr[, c(dim1_to_use, dim2_to_use, dim3_to_use)]
    dim_names <- colnames(dim_dfr)
    dim_DT <- data.table::as.data.table(dim_dfr)
    dim_DT[, cell_ID := rownames(dim_dfr)]


    ## annotated cell metadata
    cell_metadata <- combineMetadata(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        spat_enr_names = spat_enr_names
    )
    annotated_DT <- merge(cell_metadata, dim_DT, by = "cell_ID")
    spatial_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table"
    )
    if (is.null(spatial_locations)) {
        return(NULL)
    }

    annotated_DT <- merge(annotated_DT, spatial_locations, by = "cell_ID")


    if (dim_reduction_to_use == "pca") {
        pca_object <- getDimReduction(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "dimObj"
        )
        eigenvalues <- slot(pca_object, "misc")$eigenvalues

        if (!is.null(eigenvalues)) {
            total <- sum(eigenvalues)
            var_expl_vec <- (eigenvalues / total) * 100
            dim1_x_variance <- var_expl_vec[dim1_to_use]
            dim2_y_variance <- var_expl_vec[dim2_to_use]
            if (!is.null(dim3_to_use)) {
                dim3_z_variance <- var_expl_vec[3]
            }
        }
    }



    ## nn network
    if (show_NN_network) {
        # nn_network
        selected_nn_network <- getNearestNetwork(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            nn_type = nn_network_to_use,
            name = network_name,
            output = "igraph"
        )
        network_DT <- data.table::as.data.table(igraph::as_data_frame(
            selected_nn_network,
            what = "edges"
        ))

        # annotated network
        old_dim_names <- dim_names

        annotated_network_DT <- merge(network_DT, dim_DT,
            by.x = "from",
            by.y = "cell_ID"
        )
        from_dim_names <- paste0("from_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = from_dim_names
        )

        annotated_network_DT <- merge(annotated_network_DT, dim_DT,
            by.x = "to", by.y = "cell_ID"
        )
        to_dim_names <- paste0("to_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = to_dim_names
        )
    }




    ## extract spatial network
    if (show_spatial_network) {
        # Edges plus their endpoint coordinates. The network holds no
        # geometry itself, so the positions are attached here from the
        # locations -- and because that join is inner, a narrowed object
        # drops the edges whose cells are gone with no work here.
        spatial_network <- GiottoClass::annotateSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            spatial_network_name = spatial_network_name,
            spat_loc_name = spat_loc_name,
            coordinates = TRUE
        )
    } else {
        spatial_network <- NULL
    }


    ## extract spatial grid
    if (show_spatial_grid == TRUE) {
        spatial_grid <- getSpatialGrid(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = spatial_grid_name
        )
    } else {
        spatial_grid <- NULL
    }


    # create matching cell_color_code
    if (is.null(cell_color_code)) {
        if (is.character(cell_color)) {
            cell_metadata <- pDataDT(gobject,
                feat_type = feat_type,
                spat_unit = spat_unit
            )
            if (cell_color %in% colnames(cell_metadata)) {
                if (color_as_factor == TRUE) {
                    number_colors <- length(unique(cell_metadata[[cell_color]]))
                    cell_color_code <- set_default_color_discrete_cell(
                        instrs = instructions(gobject)
                    )(n = number_colors)
                    names(cell_color_code) <- unique(
                        cell_metadata[[cell_color]]
                    )
                }
            }
        }
    }


    ## subset cell selection ##
    if (!is.null(select_cells) & !is.null(select_cell_groups)) {
        if (is.null(cell_color)) {
            stop("\n selection of cells is based on cell_color paramter,
                which is a metadata column \n")
        }
        message("You have selected both individual cell IDs and a group
        of cells")
        group_cell_IDs <- annotated_DT[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- annotated_DT[get(cell_color) %in%
            select_cell_groups][["cell_ID"]]
    }


    if (!is.null(select_cells)) {
        annotated_DT_other <- annotated_DT[!annotated_DT$cell_ID %in%
            select_cells]
        annotated_DT_selected <- annotated_DT[annotated_DT$cell_ID %in%
            select_cells]

        if (show_NN_network == TRUE) {
            annotated_network_DT <- annotated_network_DT[
                annotated_network_DT$to %in% select_cells &
                    annotated_network_DT$from %in% select_cells
            ]
        }
        if (show_spatial_network == TRUE) {
            spatial_network <- spatial_network[spatial_network$to %in%
                select_cells &
                spatial_network$from %in%
                    select_cells]
        }

        # if specific cells are selected
        # annotated_DT = annotated_DT_selected
    }


    ## if no subsets are required
    if (is.null(select_cells) & is.null(select_cell_groups)) {
        annotated_DT_selected <- annotated_DT
        annotated_DT_other <- NULL
    }

    ## annotated_DT_selected = all selected cells or all cells if no selection
    ## annotated_DT_other = all not selected cells or NULL if no selection



    ########### dim plot ###########
    # 2D plot
    if (is.null(dim3_to_use)) {
        dpl <- plotly::plot_ly()
        if (show_NN_network == TRUE) {
            if (is.null(nn_network_alpha)) {
                nn_network_alpha <- 0.5
            } else if (is.character(nn_network_alpha)) {
                warning("Edge_alpha for plotly mode is not adjustable yet.
                        Default 0.5 will be set\n")
                nn_network_alpha <- 0.5
            }
            dpl <- dpl %>% plotly::add_segments(
                name = network_name,
                type = "scatter",
                x = annotated_network_DT[[from_dim_names[1]]],
                y = annotated_network_DT[[from_dim_names[2]]],
                xend = annotated_network_DT[[to_dim_names[1]]],
                yend = annotated_network_DT[[to_dim_names[2]]],
                line = list(
                    color = nn_network_color,
                    width = 0.5
                ),
                opacity = nn_network_alpha
            )
        }

        if (is.null(cell_color)) {
            # cell_color = "lightblue"
            dpl <- dpl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT_selected[[dim_names[1]]],
                y = annotated_DT_selected[[dim_names[2]]],
                # color = "lightblue",
                # colors ="lightblue",
                marker = list(
                    size = dim_point_size,
                    color = "lightblue"
                ),
                showlegend = FALSE
            )
        } else if (cell_color %in% colnames(annotated_DT_selected)) {
            if (color_as_factor) {
                annotated_DT_selected[[cell_color]] <- as.factor(
                    annotated_DT_selected[[cell_color]]
                )
            }


            dpl <- dpl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT_selected[[dim_names[1]]],
                y = annotated_DT_selected[[dim_names[2]]],
                color = annotated_DT_selected[[cell_color]],
                colors = cell_color_code,
                legendgroup = annotated_DT_selected[[cell_color]],
                marker = list(size = dim_point_size)
            )
        } else {
            stop("cell_color does not exist!\n")
        }


        if ((show_cluster_center == TRUE | show_center_label == TRUE) &
            !is.null(cell_color)) {
            annotated_DT_centers <- annotated_DT_selected[, .(
                center_1 = stats::median(get(dim_names[1])),
                center_2 = stats::median(get(dim_names[2]))
            ),
            by = cell_color
            ]
            annotated_DT_centers[[cell_color]] <- as.factor(
                annotated_DT_centers[[cell_color]]
            )
            if (show_cluster_center == TRUE) {
                dpl <- dpl %>% plotly::add_trace(
                    type = "scatter", mode = "markers",
                    x = annotated_DT_centers[["center_1"]],
                    y = annotated_DT_centers[["center_2"]],
                    color = annotated_DT_centers[[cell_color]],
                    colors = cell_color_code,
                    legendgroup = annotated_DT_centers[[cell_color]],
                    marker = list(
                        size = center_point_size, symbol = "x",
                        symbols = "x"
                    ),
                    showlegend = FALSE
                )
            }

            if (show_center_label == TRUE) {
                dpl <- dpl %>% plotly::add_text(
                    x = annotated_DT_centers[["center_1"]],
                    y = annotated_DT_centers[["center_2"]],
                    type = "scatter", mode = "text",
                    text = annotated_DT_centers[[cell_color]],
                    textposition = "middle right",
                    textfont = list(
                        color = "#000000",
                        size = label_size
                    ),
                    showlegend = FALSE
                )
            }
        }
        if (show_other_cells == TRUE) {
            dpl <- dpl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT_other[[dim_names[1]]],
                y = annotated_DT_other[[dim_names[2]]],
                marker = list(
                    size = other_point_size,
                    color = other_cell_color
                ),
                showlegend = FALSE
            )
        }
        if (dim_reduction_to_use == "pca") {
            if (!is.null(eigenvalues)) {
                x_name <- paste0("pca", "-", dim_names[1])
                y_name <- paste0("pca", "-", dim_names[2])
                x_title <- sprintf(
                    "%s explains %.02f%% of variance",
                    x_name, var_expl_vec[1]
                )
                y_title <- sprintf(
                    "%s explains %.02f%% of variance",
                    y_name, var_expl_vec[2]
                )
            }
        } else {
            x_title <- paste(dim_reduction_to_use, dim_names[1], sep = " ")
            y_title <- paste(dim_reduction_to_use, dim_names[2], sep = " ")
        }
        dpl <- dpl %>% plotly::layout(
            xaxis = list(title = x_title),
            yaxis = list(title = y_title),
            legend = list(
                x = 100, y = 0.5,
                font = list(
                    family = "sans-serif",
                    size = legend_text_size
                )
            )
        )
    }
    # 3D plot
    else if (!is.null(dim3_to_use)) {
        dpl <- plotly::plot_ly(scene = "scene1")
        if (is.null(cell_color)) {
            # cell_color = "lightblue"
            dpl <- dpl %>% plotly::add_trace(
                type = "scatter3d", mode = "markers",
                x = annotated_DT_selected[[dim_names[1]]],
                y = annotated_DT_selected[[dim_names[2]]],
                z = annotated_DT_selected[[dim_names[3]]],
                color = "lightblue",
                colors = "lightblue",
                marker = list(size = dim_point_size),
                showlegend = FALSE
            )
            # legendgroup = annotated_DT_selected[[cell_color]])
        } else {
            if (cell_color %in% colnames(annotated_DT_selected)) {
                if (is.null(cell_color_code)) {
                    number_colors <- length(unique(
                        annotated_DT_selected[[cell_color]]
                    ))
                    cell_color_code <- set_default_color_discrete_cell(
                        instrs = instructions(gobject)
                    )(n = number_colors)
                }
                if (color_as_factor) {
                    annotated_DT_selected[[cell_color]] <- as.factor(
                        annotated_DT_selected[[cell_color]]
                    )
                }
                dpl <- dpl %>% plotly::add_trace(
                    type = "scatter3d", mode = "markers",
                    x = annotated_DT_selected[[dim_names[1]]],
                    y = annotated_DT_selected[[dim_names[2]]],
                    z = annotated_DT_selected[[dim_names[3]]],
                    color = annotated_DT_selected[[cell_color]],
                    colors = cell_color_code,
                    marker = list(size = dim_point_size),
                    legendgroup = annotated_DT_selected[[cell_color]]
                )
            } else {
                stop("cell_color does not exist!\n")
            }
        }
        if (show_other_cells == TRUE) {
            dpl <- dpl %>% plotly::add_trace(
                type = "scatter3d", mode = "markers",
                x = annotated_DT_other[[dim_names[1]]],
                y = annotated_DT_other[[dim_names[2]]],
                z = annotated_DT_other[[dim_names[3]]],
                marker = list(
                    size = other_point_size,
                    color = other_cell_color
                ),
                showlegend = FALSE
            )
        }

        if (show_NN_network) {
            edges <- plotly_network(
                annotated_network_DT,
                "from_Dim.1", "from_Dim.2", "from_Dim.3",
                "to_Dim.1", "to_Dim.2", "to_Dim.3"
            )
            if (is.null(nn_network_alpha)) {
                nn_network_alpha <- 0.5
            } else if (is.character(nn_network_alpha)) {
                warning("Edge_alpha for plotly mode is not adjustable yet.
                        Default 0.5 will be set\n")
                nn_network_alpha <- 0.5
            }

            dpl <- dpl %>% plotly::add_trace(
                name = network_name,
                mode = "lines",
                type = "scatter3d",
                data = edges,
                x = ~x, y = ~y, z = ~z,
                line = list(color = nn_network_color),
                opacity = nn_network_alpha
            )
        }
        if ((show_cluster_center == TRUE | show_center_label == TRUE) &
            !is.null(cell_color)) {
            annotated_DT_centers <- annotated_DT_selected[, .(
                center_1 = stats::median(get(dim_names[1])),
                center_2 = stats::median(get(dim_names[2])),
                center_3 = stats::median(get(dim_names[3]))
            ),
            by = cell_color
            ]
            annotated_DT_centers[[cell_color]] <- as.factor(
                annotated_DT_centers[[cell_color]]
            )
            if (show_cluster_center == TRUE) {
                dpl <- dpl %>%
                    plotly::add_trace(
                        mode = "markers",
                        type = "scatter3d",
                        data = annotated_DT_centers,
                        x = ~center_1,
                        y = ~center_2,
                        z = ~center_3,
                        color = annotated_DT_centers[[cell_color]],
                        colors = cell_color_code,
                        marker = list(size = 2, symbol = "x", symbols = "x"),
                        legendgroup = annotated_DT_centers[[cell_color]],
                        showlegend = FALSE
                    )
            }
            if (show_center_label == TRUE) {
                message(" center label is not clear to see in 3D plot.
                You can shut it down with show_center_label = FALSE")
                dpl <- dpl %>%
                    plotly::add_trace(
                        mode = "text",
                        type = "scatter3d",
                        data = annotated_DT_centers,
                        x = ~center_1,
                        y = ~center_2,
                        z = ~center_3,
                        text = annotated_DT_centers[[cell_color]],
                        legendgroup = annotated_DT_centers[[cell_color]],
                        showlegend = FALSE
                    )
            }
        }
        if (dim_reduction_to_use == "pca") {
            x_name <- paste0("pca", "-", dim_names[1])
            y_name <- paste0("pca", "-", dim_names[2])
            z_name <- paste0("pca", "-", dim_names[3])
            x_title <- sprintf(
                "%s explains %.02f%% of variance",
                x_name, var_expl_vec[1]
            )
            y_title <- sprintf(
                "%s explains %.02f%% of variance",
                y_name, var_expl_vec[2]
            )
            z_title <- sprintf(
                "%s explains %.02f%% of variance",
                z_name, var_expl_vec[3]
            )
        } else {
            x_title <- paste(dim_reduction_to_use, dim_names[1], sep = " ")
            y_title <- paste(dim_reduction_to_use, dim_names[2], sep = " ")
            z_title <- paste(dim_reduction_to_use, dim_names[3], sep = " ")
        }
    }



    ############ spatial plot ##########
    if (is.null(sdimx) | is.null(sdimy)) {
        # cat('first and second dimension need to be defined, default is
        # first 2 \n')
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }

    ## 2D plot ##
    if (is.null(sdimz)) {
        spl <- plotly::plot_ly()

        if (show_spatial_network == TRUE) {
            if (is.null(spatial_network)) {
                stop("No usable spatial network specified! Please choose
                    a network with spatial_network_name=xxx")
            } else {
                if (is.null(spatial_network_alpha)) {
                    spatial_network_alpha <- 0.5
                } else if (is.character(spatial_network_alpha)) {
                    warning("Edge_alpha for plotly mode is not adjustable yet.
                            Default 0.5 will be set\n")
                    spatial_network_alpha <- 0.5
                }
                spl <- spl %>% plotly::add_segments(
                    name = spatial_network_name,
                    type = "scatter",
                    x = spatial_network[["sdimx_begin"]],
                    y = spatial_network[["sdimy_begin"]],
                    xend = spatial_network[["sdimx_end"]],
                    yend = spatial_network[["sdimy_end"]],
                    line = list(
                        color = spatial_network_color,
                        width = 0.5
                    ),
                    opacity = spatial_network_alpha
                )
            }
        }


        if (show_spatial_grid == TRUE) {
            if (is.null(spatial_grid)) {
                stop("No usable spatial grid specified! Please choose a
                    network with spatial_grid_name=xxx")
            } else {
                if (is.null(spatial_grid_color)) {
                    spatial_grid_color <- "black"
                }
                edges <- plotly_grid(spatial_grid)
                spl <- spl %>% plotly::add_segments(
                    name = "spatial_grid",
                    type = "scatter",
                    data = edges,
                    x = ~x,
                    y = ~y,
                    xend = ~x_end,
                    yend = ~y_end,
                    line = list(
                        color = spatial_grid_color,
                        width = 1
                    ),
                    opacity = spatial_grid_alpha
                )
            }
        }
        if (is.null(cell_color)) {
            # cell_color = "lightblue"
            spl <- spl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT_selected[[sdimx]],
                y = annotated_DT_selected[[sdimy]],
                # color = "lightblue",
                # colors = "lightblue",
                marker = list(
                    size = spatial_point_size,
                    color = "lightblue"
                ),
                showlegend = FALSE
            )
        } else if (cell_color %in% colnames(annotated_DT_selected)) {
            if (color_as_factor) {
                annotated_DT_selected[[cell_color]] <- as.factor(
                    annotated_DT_selected[[cell_color]]
                )
            }


            spl <- spl %>%
                plotly::add_trace(
                    type = "scatter", mode = "markers",
                    x = annotated_DT_selected[[sdimx]],
                    y = annotated_DT_selected[[sdimy]],
                    color = annotated_DT_selected[[cell_color]],
                    colors = cell_color_code,
                    legendgroup = annotated_DT_selected[[cell_color]],
                    marker = list(size = spatial_point_size),
                    showlegend = FALSE
                )
        } else {
            stop("cell_color doesn't exist!\n")
        }
        if (show_other_cells == TRUE) {
            spl <- spl %>% plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT_other[[sdimx]],
                y = annotated_DT_other[[sdimy]],
                marker = list(
                    size = other_point_size,
                    color = other_cell_color
                ),
                showlegend = FALSE
            )
        }
        spl <- spl %>% plotly::layout(
            xaxis = list(title = "X"),
            yaxis = list(title = "Y"),
            legend = list(
                x = 100, y = 0.5,
                font = list(
                    family = "sans-serif",
                    size = legend_text_size
                )
            )
        )
    }


    ## 3D plot ##
    else {
        axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))

        ratio <- plotly_axis_scale_3D(annotated_DT_selected,
            sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
            mode = axis_scale, custom_ratio = custom_ratio
        )
        spl <- plotly::plot_ly(scene = "scene2")
        if (!is.null(cell_color)) {
            if (cell_color %in% colnames(annotated_DT_selected)) {
                annotated_DT_selected[[cell_color]] <- as.factor(
                    annotated_DT_selected[[cell_color]]
                )
                spl <- spl %>%
                    plotly::add_trace(
                        type = "scatter3d", mode = "markers",
                        x = annotated_DT_selected[[sdimx]],
                        y = annotated_DT_selected[[sdimy]],
                        z = annotated_DT_selected[[sdimz]],
                        color = annotated_DT_selected[[cell_color]],
                        colors = cell_color_code,
                        legendgroup = annotated_DT_selected[[cell_color]],
                        marker = list(size = spatial_point_size),
                        showlegend = FALSE
                    )
            } else {
                stop("cell_color doesn't exist!\n")
            }
        } else {
            spl <- spl %>%
                plotly::add_trace(
                    type = "scatter3d", mode = "markers",
                    x = annotated_DT_selected$sdimx,
                    y = annotated_DT_selected$sdimy,
                    z = annotated_DT_selected$sdimz,
                    color = "lightblue",
                    colors = "lightblue",
                    # legendgroup = annotated_DT_selected[[cell_color]],
                    marker = list(size = spatial_point_size),
                    showlegend = FALSE
                )
        }
        if (show_other_cells == TRUE) {
            spl <- spl %>% plotly::add_trace(
                type = "scatter3d", mode = "markers",
                x = annotated_DT_other[[sdimx]],
                y = annotated_DT_other[[sdimy]],
                z = annotated_DT_other[[sdimz]],
                marker = list(
                    size = other_point_size,
                    color = other_cell_color
                ),
                showlegend = FALSE
            )
        }
        if (show_spatial_network == TRUE) {
            if (is.null(spatial_network)) {
                stop("No usable spatial network specified! Please choose a
                    network with spatial_network_name=xxx")
            } else {
                if (is.null(spatial_network_alpha)) {
                    spatial_network_alpha <- 0.5
                } else if (is.character(spatial_network_alpha)) {
                    warning("Edge_alpha for plotly mode is not adjustable yet.
                            Default 0.5 will be set\n")
                    spatial_network_alpha <- 0.5
                }
                edges <- plotly_network(spatial_network)

                spl <- spl %>% plotly::add_trace(
                    name = "sptial network",
                    mode = "lines",
                    type = "scatter3d",
                    data = edges,
                    x = ~x, y = ~y, z = ~z,
                    line = list(color = spatial_network_color),
                    opacity = spatial_network_alpha
                )
            }
        }

        if (show_spatial_grid == TRUE) {
            message("3D grid is not clear to view\n")
        }
    }




    if (is.null(dim3_to_use) & is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(dpl, spl,
                nrows = 2,
                titleX = TRUE, titleY = TRUE
            )
        } else {
            combo_plot <- plotly::subplot(dpl, spl,
                titleX = TRUE,
                titleY = TRUE
            )
        }
    } else if (!is.null(dim3_to_use) & is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(dpl, spl,
                nrows = 2, titleX = TRUE,
                titleY = TRUE
            ) %>%
                plotly::layout(scene = list(
                    domain = list(x = c(0, 1), y = c(0, 0.5)),
                    xaxis = list(title = x_title),
                    yaxis = list(title = y_title),
                    zaxis = list(title = z_title)
                ))
        } else {
            combo_plot <- plotly::subplot(dpl, spl,
                titleX = TRUE,
                titleY = TRUE
            ) %>%
                plotly::layout(scene = list(
                    domain = list(x = c(0, 0.5), y = c(0, 1)),
                    xaxis = list(title = x_title),
                    yaxis = list(title = y_title),
                    zaxis = list(title = z_title)
                ))
        }
    } else if (is.null(dim3_to_use) & !is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(dpl, spl,
                nrows = 2, titleX = TRUE,
                titleY = TRUE
            ) %>%
                plotly::layout(scene2 = list(
                    domain = list(x = c(0, 1), y = c(0.5, 1)),
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ))
        } else {
            combo_plot <- plotly::subplot(dpl, spl,
                titleX = TRUE,
                titleY = TRUE
            ) %>%
                plotly::layout(scene2 = list(
                    domain = list(x = c(0.5, 1), y = c(0, 1)),
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ))
        }
    } else if (!is.null(dim3_to_use) & !is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(dpl, spl,
                nrows = 2, titleX = TRUE,
                titleY = TRUE
            ) %>%
                plotly::layout(
                    scene = list(
                        domain = list(x = c(0, 1), y = c(0, 0.5)),
                        xaxis = list(title = x_title),
                        yaxis = list(title = y_title),
                        zaxis = list(title = z_title)
                    ),
                    scene2 = list(
                        domain = list(x = c(0, 1), y = c(0.5, 1)),
                        xaxis = list(title = "X", nticks = x_ticks),
                        yaxis = list(title = "Y", nticks = y_ticks),
                        zaxis = list(title = "Z", nticks = z_ticks),
                        aspectmode = "manual",
                        aspectratio = list(
                            x = ratio[[1]],
                            y = ratio[[2]],
                            z = ratio[[3]]
                        )
                    )
                )
        } else {
            combo_plot <- plotly::subplot(dpl, spl,
                titleX = TRUE,
                titleY = TRUE
            ) %>%
                plotly::layout(
                    scene = list(
                        domain = list(x = c(0, 0.5), y = c(0, 1)),
                        xaxis = list(title = x_title),
                        yaxis = list(title = y_title),
                        zaxis = list(title = z_title)
                    ),
                    scene2 = list(
                        domain = list(x = c(0.5, 1), y = c(0, 1)),
                        xaxis = list(title = "X", nticks = x_ticks),
                        yaxis = list(title = "Y", nticks = y_ticks),
                        zaxis = list(title = "Z", nticks = z_ticks),
                        aspectmode = "manual",
                        aspectratio = list(
                            x = ratio[[1]],
                            y = ratio[[2]],
                            z = ratio[[3]]
                        )
                    )
                )
        }
    }

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
        print(combo_plot)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call(
            "all_plots_save_function",
            c(list(
                gobject = gobject, plot_object = combo_plot,
                default_save_name = default_save_name
            ), save_param)
        )
    }

    ## return plot
    if (return_plot == TRUE) {
        return(combo_plot)
    }
}




# ** ####
# ** feature 3D plot ####

#' @title spatFeatPlot3D
#' @name spatFeatPlot3D
#' @description Visualize cells and gene expression according to spatial
#' coordinates
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_loc_name name of spatial locations to use
#' @param expression_values gene expression values to use
#' @param feats feats to show
#' @param spat_enr_names names of spatial enrichment results to include
#'
#' @param cluster_column cluster column to select groups
#' @param select_cell_groups select subset of cells/clusters based on
#' cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param edge_alpha alpha of edges
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#'
#' @param point_size size of point (cell)
#' @param show_legend show legend
#'
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#' @param ... additional params to pass
#' @family spatial gene expression visualizations
#' @returns plotly
#' @examples
#' g <- GiottoData::loadGiottoMini("starmap")
#' spatFeatPlot3D(g, feats = "Slc17a7")
#'
#' @export
spatFeatPlot3D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        spat_loc_name = "raw",
        expression_values = c("normalized", "scaled", "custom"),
        feats,
        spat_enr_names = NULL,
        show_network = FALSE,
        network_color = NULL,
        spatial_network_name = "Delaunay_network",
        edge_alpha = NULL,
        cluster_column = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = FALSE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        genes_high_color = NULL,
        genes_mid_color = "white",
        genes_low_color = "blue",
        show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        point_size = 2,
        show_legend = TRUE,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatFeatPlot3D",
        ...) {
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

    # data.table variables
    cell_ID <- NULL

    values <- match.arg(expression_values,
        unique(c("normalized", "scaled", "custom", expression_values))
    )

    # collect feature values using spatValues
    t_sub_expr_data_DT <- spatValues(
        gobject,
        feats = feats,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        verbose = FALSE
    )
    selected_genes <- setdiff(colnames(t_sub_expr_data_DT), "cell_ID")


    ## extract cell locations
    cell_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table"
    )
    if (is.null(cell_locations)) {
        return(NULL)
    }


    ## extract spatial network
    if (show_network) {
        # Edges plus their endpoint coordinates. The network holds no
        # geometry itself, so the positions are attached here from the
        # locations -- and because that join is inner, a narrowed object
        # drops the edges whose cells are gone with no work here.
        spatial_network <- GiottoClass::annotateSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            spatial_network_name = spatial_network_name,
            spat_loc_name = spat_loc_name,
            coordinates = TRUE
        )
    } else {
        spatial_network <- NULL
    }

    ## extract spatial grid
    if (show_grid) {
        spatial_grid <- getSpatialGrid(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = spatial_grid_name
        )
    } else {
        spatial_grid <- NULL
    }

    ## extract cell metadata
    cell_metadata <- try(
        expr = combineMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            spat_loc_name = spat_loc_name,
            spat_enr_names = spat_enr_names
        ),
        silent = TRUE
    )


    if (inherits(cell_metadata, "try-error")) {
        cell_locations_metadata <- cell_locations
    } else if (nrow(cell_metadata) == 0) {
        cell_locations_metadata <- cell_locations
    } else {
        cell_locations_metadata <- cell_metadata
    }


    if (!is.null(select_cells) & !is.null(select_cell_groups)) {
        message("You have selected both individual cell IDs and a group
        of cells")
        group_cell_IDs <- cell_locations_metadata[get(cluster_column) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- cell_locations_metadata[get(cluster_column) %in%
            select_cell_groups][["cell_ID"]]
    }

    if (!is.null(select_cells)) {
        cell_locations_metadata_other <- cell_locations_metadata[
            !cell_locations_metadata$cell_ID %in% select_cells
        ]
        cell_locations_metadata_selected <- cell_locations_metadata[
            cell_locations_metadata$cell_ID %in% select_cells
        ]
        spatial_network <- spatial_network[spatial_network$to %in%
            select_cells & spatial_network$from %in% select_cells]

        # if specific cells are selected
        cell_locations_metadata <- cell_locations_metadata_selected
    }

    cell_locations_metadata_genes <- merge(cell_locations_metadata,
        t_sub_expr_data_DT,
        by = "cell_ID"
    )



    ## plotting ##
    axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))

    ratio <- plotly_axis_scale_3D(cell_locations_metadata_genes,
        sdimx = "sdimx", sdimy = "sdimy", sdimz = "sdimz",
        mode = axis_scale, custom_ratio = custom_ratio
    )


    ## spatial network data
    if (!is.null(spatial_network) & show_network == TRUE) {
        edges <- plotly_network(spatial_network)
    }

    ## Point layer
    if (length(selected_genes) > 4) {
        stop("\n The max number of genes showed together is 4.Otherwise
            it will be too small to see\n
            \n If you have more genes to show, please divide them
            into groups\n")
    }
    savelist <- list()
    for (i in seq_len(length(selected_genes))) {
        gene <- selected_genes[i]
        if (!is.null(genes_high_color)) {
            if (length(genes_high_color) != length(selected_genes) &
                length(genes_high_color) != 1) {
                stop("\n The number of genes and their corresbonding do
                    not match\n")
            } else if (length(genes_high_color) == 1) {
                genes_high_color <- rep(
                    genes_high_color,
                    length(selected_genes)
                )
            }
        } else {
            genes_high_color <- rep("red", length(selected_genes))
        }
        pl <- plotly::plot_ly(
            name = gene,
            scene = paste("scene", i, sep = "")
        ) %>%
            plotly::add_trace(
                data = cell_locations_metadata_genes,
                type = "scatter3d", mode = "markers",
                x = ~sdimx, y = ~sdimy, z = ~sdimz,
                marker = list(size = point_size),
                color = cell_locations_metadata_genes[[gene]],
                colors = c(
                    genes_low_color, genes_mid_color,
                    genes_high_color[i]
                )
            )

        if (show_other_cells == TRUE) {
            pl <- pl %>% plotly::add_trace(
                name = "unselected cells",
                data = cell_locations_metadata_other,
                type = "scatter3d", mode = "markers",
                x = ~sdimx, y = ~sdimy, z = ~sdimz,
                marker = list(size = other_point_size, color = other_cell_color)
            )
        }


        ## plot spatial network
        if (show_network == TRUE) {
            if (is.null(network_color)) {
                network_color <- "lightblue"
            }
            if (is.null(edge_alpha)) {
                edge_alpha <- 0.5
            } else if (is.character(edge_alpha)) {
                edge_alpha <- 0.5
                message("Edge_alpha for plotly mode is not adjustable yet.
                Default 0.5 will be set")
            }
            pl <- pl %>% plotly::add_trace(
                name = "sptial network",
                mode = "lines",
                type = "scatter3d",
                data = edges,
                x = ~x, y = ~y, z = ~z,
                line = list(color = network_color, width = 0.5),
                opacity = edge_alpha,
                showlegend = FALSE
            )
        }


        ## plot spatial grid
        if (!is.null(spatial_grid) & show_grid == TRUE) {
            message("spatial grid is not clear in 3D plot")
        }

        pl <- pl %>% plotly::colorbar(title = gene)
        savelist[[gene]] <- pl
    }


    if (length(savelist) == 1) {
        cowplot <- savelist[[1]] %>% plotly::layout(scene = list(
            xaxis = list(title = "X", nticks = x_ticks),
            yaxis = list(title = "Y", nticks = y_ticks),
            zaxis = list(title = "Z", nticks = z_ticks),
            aspectmode = "manual",
            aspectratio = list(
                x = ratio[[1]],
                y = ratio[[2]],
                z = ratio[[3]]
            )
        ))
    } else if (length(savelist) == 2) {
        cowplot <- plotly::subplot(savelist) %>%
            plotly::layout(
                scene = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                scene2 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                # annotations = annotations,
                legend = list(x = 100, y = 0)
            )
    } else if (length(savelist) == 3) {
        cowplot <- plotly::subplot(savelist) %>%
            plotly::layout(
                scene = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                scene2 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                scene3 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                legend = list(x = 100, y = 0)
            )
    } else if (length(savelist) == 4) {
        cowplot <- plotly::subplot(savelist) %>%
            plotly::layout(
                scene = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                scene2 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                scene3 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                scene4 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ),
                legend = list(x = 100, y = 0)
            )
    }


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
        print(cowplot)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call(
            "all_plots_save_function",
            c(list(
                gobject = gobject, plot_object = cowplot,
                default_save_name = default_save_name
            ), save_param)
        )
    }

    ## return plot
    if (return_plot == TRUE) {
        return(cowplot)
    }
}

#' @describeIn spatFeatPlot3D deprecated
#' @export
spatGenePlot3D <- function(...) {
    deprecate_warn(
        when = "0.2.0",
        what = "spatGenePlot3D()",
        with = "spatFeatPlot3D()"
    )
    spatFeatPlot3D(...)
}


#' @title dimFeatPlot3D
#' @name dimFeatPlot3D
#' @description Visualize cells and gene expression according to
#' dimension reduction coordinates
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#'
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param network_color color of NN network
#'
#' @param cluster_column cluster column to select groups
#' @param select_cell_groups select subset of cells/clusters based on
#' cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#'
#' @param genes_high_color color for high expression levels
#' @param genes_mid_color color for medium expression levels
#' @param genes_low_color color for low expression levels
#'
#' @param show_legend show legend
#' @details Description of parameters.
#' @family dimension reduction gene expression visualizations
#' @returns plotly
#' @examples
#' g <- GiottoData::loadGiottoMini("starmap")
#'
#' dimFeatPlot3D(g, genes = "Slc17a7", dim_reduction_name = "3D_umap")
#' @export
dimFeatPlot3D <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        genes = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        dim3_to_use = 3,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        network_color = "lightgray",
        cluster_column = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        edge_alpha = NULL,
        point_size = 2,
        genes_high_color = NULL,
        genes_mid_color = "white",
        genes_low_color = "blue",
        show_legend = TRUE,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "dimFeatPlot3D") {
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

    # data.table variables
    cell_ID <- NULL

    ## select genes ##
    values <- match.arg(expression_values, c("normalized", "scaled", "custom"))

    # collect feature values using spatValues
    t_sub_expr_data_DT <- spatValues(
        gobject,
        feats = genes,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        verbose = FALSE
    )
    selected_genes <- setdiff(colnames(t_sub_expr_data_DT), "cell_ID")


    ## dimension reduction ##
    dim_dfr <- getDimReduction(gobject,
        reduction = "cells",
        reduction_method = dim_reduction_to_use,
        name = dim_reduction_name,
        output = "data.table"
    )
    dim_dfr <- dim_dfr[, c(dim1_to_use, dim2_to_use, dim3_to_use)]
    dim_names <- colnames(dim_dfr)
    dim_DT <- data.table::as.data.table(dim_dfr)
    dim_DT[, cell_ID := rownames(dim_dfr)]

    ## annotated cell metadata
    cell_metadata <- pDataDT(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit
    )
    annotated_DT <- merge(cell_metadata, dim_DT, by = "cell_ID")



    # create input for network
    if (show_NN_network == TRUE) {
        # nn_network
        selected_nn_network <- getNearestNetwork(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            nn_type = nn_network_to_use,
            name = network_name,
            output = "igraph"
        )
        network_DT <- data.table::as.data.table(igraph::as_data_frame(
            selected_nn_network,
            what = "edges"
        ))

        # annotated network
        old_dim_names <- dim_names

        annotated_network_DT <- merge(
            network_DT, dim_DT,
            by.x = "from", by.y = "cell_ID"
        )
        from_dim_names <- paste0("from_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = from_dim_names
        )

        annotated_network_DT <- merge(
            annotated_network_DT, dim_DT,
            by.x = "to", by.y = "cell_ID"
        )
        to_dim_names <- paste0("to_", old_dim_names)
        data.table::setnames(annotated_network_DT,
            old = old_dim_names,
            new = to_dim_names
        )
    }


    ## create subsets if needed
    if (!is.null(select_cells) & !is.null(select_cell_groups)) {
        if (is.null(cluster_column)) {
            stop("\n selection of cells is based on cell_color paramter,
                which is a metadata column \n")
        }
        message("You have selected both individual cell IDs and a group
        of cells")
        group_cell_IDs <- annotated_DT[get(cluster_column) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- annotated_DT[get(cluster_column) %in%
            select_cell_groups][["cell_ID"]]
    }

    if (!is.null(select_cells)) {
        annotated_DT_other <- annotated_DT[
            !annotated_DT$cell_ID %in% select_cells
        ]
        annotated_DT_selected <- annotated_DT[
            annotated_DT$cell_ID %in% select_cells
        ]

        if (show_NN_network == TRUE) {
            annotated_network_DT <- annotated_network_DT[
                annotated_network_DT$to %in% select_cells &
                    annotated_network_DT$from %in% select_cells
            ]
        }

        # if specific cells are selected
        annotated_DT <- annotated_DT_selected
    }

    ## if no subsets are required
    if (is.null(select_cells) & is.null(select_cell_groups)) {
        annotated_DT_selected <- annotated_DT
        annotated_DT_other <- NULL
    }

    ## merge gene info
    annotated_gene_DT <- merge(annotated_DT, t_sub_expr_data_DT, by = "cell_ID")



    ## visualize multipe plots ##
    ## 3D plots ##


    if (show_NN_network == TRUE) {
        edges <- plotly_network(
            annotated_network_DT,
            "from_Dim.1", "from_Dim.2", "from_Dim.3",
            "to_Dim.1", "to_Dim.2", "to_Dim.3"
        )
    }
    ## Point layer
    if (length(selected_genes) > 4) {
        stop("\n The max number of genes showed together is 4.Otherwise
            it will be too small to see\n
            \n If you have more genes to show, please divide them into
            groups\n")
    }
    if (!is.null(genes_high_color)) {
        if (length(genes_high_color) != length(selected_genes) &
            length(genes_high_color) != 1) {
            stop("\n The number of genes and their corresbonding do not
                match\n")
        }
    } else if (is.null(genes_high_color)) {
        genes_high_color <- rep("red", length(selected_genes))
    } else {
        genes_high_color <- rep(genes_high_color, length(selected_genes))
    }

    titleX <- title <- paste(dim_reduction_to_use, dim_names[1], sep = " ")
    titleY <- title <- paste(dim_reduction_to_use, dim_names[2], sep = " ")
    titleZ <- title <- paste(dim_reduction_to_use, dim_names[3], sep = " ")
    savelist <- list()
    for (i in seq_len(length(selected_genes))) {
        gene <- selected_genes[i]

        pl <- plotly::plot_ly(name = gene, scene = paste("scene", i, sep = ""))
        pl <- pl %>% plotly::add_trace(
            data = annotated_gene_DT, type = "scatter3d", mode = "markers",
            x = annotated_gene_DT[[dim_names[1]]],
            y = annotated_gene_DT[[dim_names[2]]],
            z = annotated_gene_DT[[dim_names[3]]],
            color = annotated_gene_DT[[gene]],
            colors = c(genes_low_color, genes_mid_color, genes_high_color[i]),
            marker = list(size = point_size)
        )
        if (show_other_cells == TRUE) {
            pl <- pl %>% plotly::add_trace(
                name = "unselected cells",
                data = annotated_DT_other,
                type = "scatter3d", mode = "markers",
                x = annotated_DT_other[[dim_names[1]]],
                y = annotated_DT_other[[dim_names[2]]],
                z = annotated_DT_other[[dim_names[3]]],
                marker = list(size = other_point_size, color = other_cell_color)
            )
        }

        ## plot spatial network
        if (show_NN_network == TRUE) {
            pl <- pl %>% plotly::add_trace(
                name = "sptial network", mode = "lines",
                type = "scatter3d", opacity = edge_alpha,
                showlegend = FALSE,
                data = edges,
                x = ~x, y = ~y, z = ~z,
                line = list(
                    color = network_color,
                    width = 0.5
                )
            )
        }
        pl <- pl %>% plotly::colorbar(title = gene)
        savelist[[gene]] <- pl
    }

    if (length(savelist) == 1) {
        cowplot <- savelist[[1]] %>% plotly::layout(scene = list(
            xaxis = list(title = titleX),
            yaxis = list(title = titleY),
            zaxis = list(title = titleZ)
        ))
    } else if (length(savelist) == 2) {
        cowplot <- plotly::subplot(
            savelist,
            titleX = TRUE, titleY = TRUE
        ) %>%
            plotly::layout(
                scene = list(
                    domain = list(x = c(0, 0.5), y = c(0, 1)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                scene2 = list(
                    domain = list(x = c(0.5, 1), y = c(0, 1)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                legend = list(x = 100, y = 0)
            )
    } else if (length(savelist) == 3) {
        cowplot <- plotly::subplot(
            savelist,
            titleX = TRUE, titleY = TRUE
        ) %>%
            plotly::layout(
                scene = list(
                    domain = list(x = c(0, 0.5), y = c(0, 0.5)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                scene2 = list(
                    domain = list(x = c(0.5, 1), y = c(0, 0.5)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                scene3 = list(
                    domain = list(x = c(0, 0.5), y = c(0.5, 1)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                legend = list(x = 100, y = 0)
            )
    } else if (length(savelist) == 4) {
        cowplot <- plotly::subplot(savelist) %>%
            plotly::layout(
                scene = list(
                    domain = list(x = c(0, 0.5), y = c(0, 0.5)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                scene2 = list(
                    domain = list(x = c(0.5, 1), y = c(0, 0.5)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                scene3 = list(
                    domain = list(x = c(0, 0.5), y = c(0.5, 1)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                scene4 = list(
                    domain = list(x = c(0.5, 1), y = c(0.5, 1)),
                    xaxis = list(title = titleX),
                    yaxis = list(title = titleY),
                    zaxis = list(title = titleZ)
                ),
                legend = list(x = 100, y = 0)
            )
    }

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
        print(cowplot)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call(
            "all_plots_save_function",
            c(list(
                gobject = gobject, plot_object = cowplot,
                default_save_name = default_save_name
            ), save_param)
        )
    }

    ## return plot
    if (return_plot == TRUE) {
        return(cowplot)
    }
}

#' @describeIn dimFeatPlot3D deprecated
#' @param ... additional params to pass
#' @export
dimGenePlot3D <- function(...) {
    deprecate_warn(
        when = "0.2.0",
        what = "dimGenePlot3D()",
        with = "dimFeatPlot3D()"
    )
    dimFeatPlot3D(...)
}



#' @title spatDimFeatPlot3D
#' @name spatDimFeatPlot3D
#' @description Visualize cells according to spatial AND dimension
#' reduction coordinates in ggplot mode
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_loc_name name of spatial locations to use
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param sdimx spatial dimension to use on x-axis
#' @param sdimy spatial dimension to use on y-axis
#' @param sdimz spatial dimension to use on z-axis
#' @param genes genes to show
#'
#' @param cluster_column cluster column to select groups
#' @param select_cell_groups select subset of cells/clusters based on
#' cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param dim_point_size dim reduction plot: point size
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param nn_network_color color of NN network
#' @param nn_network_alpha alpha of NN network
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#'
#' @param label_size size of labels
#' @param genes_high_color color for high expression levels
#' @param genes_mid_color color for medium expression levels
#' @param genes_low_color color for low expression levels
#'
#' @param show_spatial_network show spatial network (boolean)
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param spatial_network_alpha alpha of spatial network
#'
#' @param show_spatial_grid show spatial grid (boolean)
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_grid_alpha alpha of spatial grid
#'
#' @param spatial_point_size spatial plot: point size
#' @param legend_text_size size of legend
#'
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#' @details Description of parameters.
#' @family spatial and dimension reduction gene expression visualizations
#' @returns plotly
#' @examples
#' g <- GiottoData::loadGiottoMini("starmap")
#' spatDimFeatPlot3D(g, genes = "Slc17a7")
#'
#' @export
spatDimFeatPlot3D <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        spat_loc_name = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        plot_alignment = c("horizontal", "vertical"),
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        dim3_to_use = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        sdimz = "sdimz",
        genes,
        cluster_column = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1.5,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        nn_network_color = "lightgrey",
        nn_network_alpha = 0.5,
        network_name = "sNN.pca",
        label_size = 16,
        genes_low_color = "blue",
        genes_mid_color = "white",
        genes_high_color = "red",
        dim_point_size = 3,
        show_spatial_network = FALSE,
        spatial_network_name = "Delaunay_network",
        spatial_network_color = "lightgray",
        spatial_network_alpha = 0.5,
        show_spatial_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        spatial_grid_color = NULL,
        spatial_grid_alpha = 0.5,
        spatial_point_size = 3,
        legend_text_size = 12,
        axis_scale = c("cube", "real", "custom"),
        custom_ratio = NULL,
        x_ticks = NULL,
        y_ticks = NULL,
        z_ticks = NULL,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatDimFeatPlot3D") {
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

    # data.table variables
    cell_ID <- NULL

    plot_alignment <- match.arg(plot_alignment,
        choices = c("horizontal", "vertical")
    )

    ########### data prepare ###########
    ## select genes ##
    if (length(genes) > 1) {
        warning("\n Now 3D mode can just accept one gene, only the first
                gene will be plot\n")
        genes <- genes[1]
    }
    values <- match.arg(expression_values, c("normalized", "scaled", "custom"))

    # collect feature values using spatValues
    t_sub_expr_data_DT <- spatValues(
        gobject,
        feats = genes,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        verbose = FALSE
    )
    selected_genes <- setdiff(colnames(t_sub_expr_data_DT), "cell_ID")


    ## dimension reduction ##
    dim_dfr <- getDimReduction(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        reduction = "cells",
        reduction_method = dim_reduction_to_use,
        name = dim_reduction_name,
        output = "matrix"
    )
    dim_dfr <- dim_dfr[, c(dim1_to_use, dim2_to_use, dim3_to_use)]
    dim_names <- colnames(dim_dfr)
    dim_DT <- data.table::as.data.table(dim_dfr)
    dim_DT[, cell_ID := rownames(dim_dfr)]


    ## annotated cell metadata
    cell_metadata <- pDataDT(gobject,
        feat_type = feat_type,
        spat_unit = spat_unit
    )
    cell_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table"
    )
    if (is.null(cell_locations)) {
        return(NULL)
    }

    annotated_DT <- merge(cell_metadata, dim_DT, by = "cell_ID")
    annotated_DT <- merge(annotated_DT, cell_locations, by = "cell_ID")
    annotated_DT <- merge(annotated_DT, t_sub_expr_data_DT, by = "cell_ID")


    ## nn network
    if (show_NN_network) {
        # nn_network
        selected_nn_network <- getNearestNetwork(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            nn_type = nn_network_to_use,
            name = network_name,
            output = "igraph"
        )
        network_DT <- data.table::as.data.table(igraph::as_data_frame(
            selected_nn_network,
            what = "edges"
        ))

        # annotated network
        old_dim_names <- dim_names

        annotated_network_DT <- merge(
            network_DT, dim_DT,
            by.x = "from", by.y = "cell_ID"
        )
        from_dim_names <- paste0("from_", old_dim_names)
        data.table::setnames(
            annotated_network_DT,
            old = old_dim_names, new = from_dim_names
        )

        annotated_network_DT <- merge(annotated_network_DT, dim_DT,
            by.x = "to", by.y = "cell_ID"
        )
        to_dim_names <- paste0("to_", old_dim_names)
        data.table::setnames(
            annotated_network_DT,
            old = old_dim_names, new = to_dim_names
        )
    }


    ## extract spatial network
    if (show_spatial_network == TRUE) {
        # Edges plus their endpoint coordinates. The network holds no
        # geometry itself, so the positions are attached here from the
        # locations -- and because that join is inner, a narrowed object
        # drops the edges whose cells are gone with no work here.
        spatial_network <- GiottoClass::annotateSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            spatial_network_name = spatial_network_name,
            spat_loc_name = spat_loc_name,
            coordinates = TRUE
        )
    } else {
        spatial_network <- NULL
    }

    ## extract spatial grid
    if (show_spatial_grid == TRUE) {
        spatial_grid <- getSpatialGrid(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            name = spatial_grid_name
        )
    } else {
        spatial_grid <- NULL
    }


    ## select subset of cells ##
    if (!is.null(select_cells) & !is.null(select_cell_groups)) {
        if (is.null(cluster_column)) {
            stop("\n selection of cells is based on cell_color paramter,
                which is a metadata column \n")
        }
        message("You have selected both individual cell IDs and a group
        of cells")
        group_cell_IDs <- annotated_DT[get(cluster_column) %in%
            select_cell_groups][["cell_ID"]]
        select_cells <- unique(c(select_cells, group_cell_IDs))
    } else if (!is.null(select_cell_groups)) {
        select_cells <- annotated_DT[get(cluster_column) %in%
            select_cell_groups][["cell_ID"]]
    }

    if (!is.null(select_cells)) {
        annotated_DT_other <- annotated_DT[
            !annotated_DT$cell_ID %in% select_cells
        ]
        annotated_DT_selected <- annotated_DT[
            annotated_DT$cell_ID %in% select_cells
        ]

        if (show_NN_network == TRUE) {
            annotated_network_DT <- annotated_network_DT[
                annotated_network_DT$to %in% select_cells &
                    annotated_network_DT$from %in% select_cells
            ]
        }
        if (show_spatial_network == TRUE) {
            spatial_network <- spatial_network[
                spatial_network$to %in% select_cells &
                    spatial_network$from %in% select_cells
            ]
        }

        # if specific cells are selected
        annotated_DT <- annotated_DT_selected
    }

    ## if no subsets are required
    if (is.null(select_cells) & is.null(select_cell_groups)) {
        annotated_DT_selected <- annotated_DT
        annotated_DT_other <- NULL
    }




    ########### dim plot ###########
    # 2D plot
    if (is.null(dim3_to_use)) {
        dpl <- plotly::plot_ly()
        if (show_NN_network == TRUE) {
            if (is.null(nn_network_alpha)) {
                nn_network_alpha <- 0.5
            } else if (is.character(nn_network_alpha)) {
                warning("Edge_alpha for plotly mode is not adjustable yet.
                        Default 0.5 will be set\n")
                nn_network_alpha <- 0.5
            }
            dpl <- dpl %>% plotly::add_segments(
                name = network_name,
                type = "scatter",
                x = annotated_network_DT[[from_dim_names[1]]],
                y = annotated_network_DT[[from_dim_names[2]]],
                xend = annotated_network_DT[[to_dim_names[1]]],
                yend = annotated_network_DT[[to_dim_names[2]]],
                line = list(
                    color = nn_network_color,
                    width = 0.5
                ),
                opacity = nn_network_alpha
            )
        }

        dpl <- dpl %>%
            plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT[[dim_names[1]]],
                y = annotated_DT[[dim_names[2]]],
                color = annotated_DT[[selected_genes]],
                colors = c(
                    genes_low_color, genes_mid_color,
                    genes_high_color
                ),
                marker = list(size = dim_point_size),
                showlegend = FALSE
            )

        if (show_other_cells == TRUE) {
            dpl <- dpl %>%
                plotly::add_trace(
                    type = "scatter", mode = "markers",
                    x = annotated_DT_other[[dim_names[1]]],
                    y = annotated_DT_other[[dim_names[2]]],
                    marker = list(
                        size = other_point_size,
                        color = other_cell_color
                    ),
                    showlegend = FALSE
                )
        }

        x_title <- paste(dim_reduction_to_use, dim_names[1], sep = " ")
        y_title <- paste(dim_reduction_to_use, dim_names[2], sep = " ")

        dpl <- dpl %>% plotly::layout(
            xaxis = list(title = x_title),
            yaxis = list(title = y_title),
            legend = list(
                x = 100, y = 0.5,
                font = list(
                    family = "sans-serif",
                    size = legend_text_size
                )
            )
        )
    }
    # 3D plot
    else if (!is.null(dim3_to_use)) {
        dpl <- plotly::plot_ly(scene = "scene1")

        dpl <- dpl %>% plotly::add_trace(
            type = "scatter3d", mode = "markers",
            x = annotated_DT[[dim_names[1]]],
            y = annotated_DT[[dim_names[2]]],
            z = annotated_DT[[dim_names[3]]],
            color = annotated_DT[[selected_genes]],
            colors = c(genes_low_color, genes_mid_color, genes_high_color),
            marker = list(size = dim_point_size),
            showlegend = FALSE
        )
        # legendgroup = annotated_DT[[cell_color]])
        if (show_other_cells == TRUE) {
            dpl <- dpl %>% plotly::add_trace(
                type = "scatter3d", mode = "markers",
                x = annotated_DT_other[[dim_names[1]]],
                y = annotated_DT_other[[dim_names[2]]],
                z = annotated_DT_other[[dim_names[3]]],
                marker = list(
                    size = other_point_size,
                    color = other_cell_color
                ),
                showlegend = FALSE
            )
        }

        if (show_NN_network) {
            edges <- plotly_network(
                annotated_network_DT,
                "from_Dim.1", "from_Dim.2", "from_Dim.3",
                "to_Dim.1", "to_Dim.2", "to_Dim.3"
            )
            if (is.null(nn_network_alpha)) {
                nn_network_alpha <- 0.5
            } else if (is.character(nn_network_alpha)) {
                warning("Edge_alpha for plotly mode is not adjustable yet.
                        Default 0.5 will be set\n")
                nn_network_alpha <- 0.5
            }

            dpl <- dpl %>% plotly::add_trace(
                name = network_name,
                mode = "lines",
                type = "scatter3d",
                data = edges,
                x = ~x, y = ~y, z = ~z,
                line = list(color = nn_network_color),
                opacity = nn_network_alpha
            )
        }


        x_title <- paste(dim_reduction_to_use, dim_names[1], sep = " ")
        y_title <- paste(dim_reduction_to_use, dim_names[2], sep = " ")
        z_title <- paste(dim_reduction_to_use, dim_names[3], sep = " ")
    }
    dpl <- dpl %>% plotly::colorbar(title = selected_genes)


    ########### spatial plot ###########
    if (is.null(sdimx) | is.null(sdimy)) {
        # cat('first and second dimenion need to be defined,
        # default is first 2 \n')
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }

    # 2D plot
    if (is.null(sdimz)) {
        spl <- plotly::plot_ly()

        if (show_spatial_network == TRUE) {
            if (is.null(spatial_network)) {
                stop("No usable spatial network specified! Please choose a
                    network with spatial_network_name=xxx")
            } else {
                if (is.null(spatial_network_alpha)) {
                    spatial_network_alpha <- 0.5
                } else if (is.character(spatial_network_alpha)) {
                    warning("Edge_alpha for plotly mode is not adjustable yet.
                            Default 0.5 will be set\n")
                    spatial_network_alpha <- 0.5
                }
                spl <- spl %>% plotly::add_segments(
                    name = spatial_network_name,
                    type = "scatter",
                    x = spatial_network[["sdimx_begin"]],
                    y = spatial_network[["sdimy_begin"]],
                    xend = spatial_network[["sdimx_end"]],
                    yend = spatial_network[["sdimy_end"]],
                    line = list(
                        color = spatial_network_color,
                        width = 0.5
                    ),
                    opacity = spatial_network_alpha
                )
            }
        }
        if (show_spatial_grid == TRUE) {
            if (is.null(spatial_grid)) {
                stop("No usable spatial grid specified! Please choose a
                    network with spatial_grid_name=xxx")
            } else {
                if (is.null(spatial_grid_color)) {
                    spatial_grid_color <- "black"
                }
                edges <- plotly_grid(spatial_grid)
                spl <- spl %>% plotly::add_segments(
                    name = "spatial_grid",
                    type = "scatter",
                    data = edges,
                    x = ~x,
                    y = ~y,
                    xend = ~x_end,
                    yend = ~y_end,
                    line = list(
                        color = spatial_grid_color,
                        width = 1
                    ),
                    opacity = spatial_grid_alpha
                )
            }
        }

        spl <- spl %>%
            plotly::add_trace(
                type = "scatter", mode = "markers",
                x = annotated_DT[[sdimx]],
                y = annotated_DT[[sdimy]],
                color = annotated_DT[[selected_genes]],
                colors = c(
                    genes_low_color, genes_mid_color,
                    genes_high_color
                ),
                marker = list(size = spatial_point_size),
                showlegend = FALSE
            )
        if (show_other_cells == TRUE) {
            spl <- spl %>%
                plotly::add_trace(
                    type = "scatter", mode = "markers",
                    x = annotated_DT_other[[sdimx]],
                    y = annotated_DT_other[[sdimy]],
                    marker = list(
                        size = other_point_size,
                        color = other_cell_color
                    ),
                    showlegend = FALSE
                )
        }

        spl <- spl %>% plotly::layout(
            xaxis = list(title = "X"),
            yaxis = list(title = "Y"),
            legend = list(
                x = 100, y = 0.5,
                font = list(
                    family = "sans-serif",
                    size = legend_text_size
                )
            )
        )
    }


    # 3D plot
    else {
        axis_scale <- match.arg(axis_scale, c("cube", "real", "custom"))
        ratio <- plotly_axis_scale_3D(annotated_DT,
            sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
            mode = axis_scale, custom_ratio = custom_ratio
        )


        spl <- plotly::plot_ly(scene = "scene2")

        spl <- spl %>%
            plotly::add_trace(
                type = "scatter3d", mode = "markers",
                x = annotated_DT[[sdimx]],
                y = annotated_DT[[sdimy]],
                z = annotated_DT[[sdimz]],
                color = annotated_DT[[selected_genes]],
                colors = c(
                    genes_low_color, genes_mid_color,
                    genes_high_color
                ),
                # legendgroup = annotated_DT[[cell_color]],
                marker = list(size = spatial_point_size),
                showlegend = FALSE
            )
        if (show_other_cells == TRUE) {
            spl <- spl %>%
                plotly::add_trace(
                    type = "scatter3d", mode = "markers",
                    x = annotated_DT_other[[sdimx]],
                    y = annotated_DT_other[[sdimy]],
                    z = annotated_DT_other[[sdimz]],
                    marker = list(
                        size = other_point_size,
                        color = other_cell_color
                    ),
                    showlegend = FALSE
                )
        }

        if (show_spatial_network == TRUE) {
            if (is.null(spatial_network)) {
                stop("No usable spatial network specified! Please choose a
                    network with spatial_network_name=xxx")
            } else {
                if (is.null(spatial_network_alpha)) {
                    spatial_network_alpha <- 0.5
                } else if (is.character(spatial_network_alpha)) {
                    warning("Edge_alpha for plotly mode is not adjustable yet.
                            Default 0.5 will be set\n")
                    spatial_network_alpha <- 0.5
                }
                edges <- plotly_network(spatial_network)

                spl <- spl %>% plotly::add_trace(
                    name = "sptial network",
                    mode = "lines",
                    type = "scatter3d",
                    data = edges,
                    x = ~x, y = ~y, z = ~z,
                    line = list(color = spatial_network_color),
                    opacity = spatial_network_alpha
                )
            }
        }

        if (show_spatial_grid == TRUE) {
            message("3D grid is not clear to view")
        }
    }



    spl <- plotly::hide_colorbar(spl)
    if (is.null(dim3_to_use) & is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(
                dpl, spl,
                nrows = 2, titleX = TRUE, titleY = TRUE
            )
        } else {
            combo_plot <- plotly::subplot(
                dpl, spl,
                titleX = TRUE, titleY = TRUE
            )
        }
    } else if (!is.null(dim3_to_use) & is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(
                dpl, spl,
                nrows = 2, titleX = TRUE, titleY = TRUE
            ) %>%
                plotly::layout(scene = list(
                    domain = list(x = c(0, 1), y = c(0, 0.5)),
                    xaxis = list(title = x_title),
                    yaxis = list(title = y_title),
                    zaxis = list(title = z_title)
                ))
        } else {
            combo_plot <- plotly::subplot(
                dpl, spl,
                titleX = TRUE, titleY = TRUE
            ) %>%
                plotly::layout(scene = list(
                    domain = list(x = c(0, 0.5), y = c(0, 1)),
                    xaxis = list(title = x_title),
                    yaxis = list(title = y_title),
                    zaxis = list(title = z_title)
                ))
        }
    } else if (is.null(dim3_to_use) & !is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(
                dpl, spl,
                nrows = 2, titleX = TRUE, titleY = TRUE
            ) %>%
                plotly::layout(scene2 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ))
        } else {
            combo_plot <- plotly::subplot(
                dpl, spl,
                titleX = TRUE, titleY = TRUE
            ) %>%
                plotly::layout(scene2 = list(
                    xaxis = list(title = "X", nticks = x_ticks),
                    yaxis = list(title = "Y", nticks = y_ticks),
                    zaxis = list(title = "Z", nticks = z_ticks),
                    aspectmode = "manual",
                    aspectratio = list(
                        x = ratio[[1]],
                        y = ratio[[2]],
                        z = ratio[[3]]
                    )
                ))
        }
    } else if (!is.null(dim3_to_use) & !is.null(sdimz)) {
        if (plot_alignment == "vertical") {
            combo_plot <- plotly::subplot(
                dpl, spl,
                nrows = 2, titleX = TRUE, titleY = TRUE
            ) %>%
                plotly::layout(
                    scene = list(
                        domain = list(x = c(0, 1), y = c(0, 0.5)),
                        xaxis = list(title = x_title),
                        yaxis = list(title = y_title),
                        zaxis = list(title = z_title)
                    ),
                    scene2 = list(
                        xaxis = list(title = "X", nticks = x_ticks),
                        yaxis = list(title = "Y", nticks = y_ticks),
                        zaxis = list(title = "Z", nticks = z_ticks),
                        aspectmode = "manual",
                        aspectratio = list(
                            x = ratio[[1]],
                            y = ratio[[2]],
                            z = ratio[[3]]
                        )
                    )
                )
        } else {
            combo_plot <- plotly::subplot(
                dpl, spl,
                titleX = TRUE, titleY = TRUE
            ) %>%
                plotly::layout(
                    scene = list(
                        domain = list(x = c(0, 0.5), y = c(0, 1)),
                        xaxis = list(title = x_title),
                        yaxis = list(title = y_title),
                        zaxis = list(title = z_title)
                    ),
                    scene2 = list(
                        xaxis = list(title = "X", nticks = x_ticks),
                        yaxis = list(title = "Y", nticks = y_ticks),
                        zaxis = list(title = "Z", nticks = z_ticks),
                        aspectmode = "manual",
                        aspectratio = list(
                            x = ratio[[1]],
                            y = ratio[[2]],
                            z = ratio[[3]]
                        )
                    )
                )
        }
    }

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
        print(combo_plot)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call(
            "all_plots_save_function",
            c(list(
                gobject = gobject, plot_object = combo_plot,
                default_save_name = default_save_name
            ), save_param)
        )
    }

    ## return plot
    if (return_plot == TRUE) {
        return(combo_plot)
    }
}

#' @describeIn spatDimFeatPlot3D deprecated
#' @param ... additional params to pass
#' @export
spatDimGenePlot3D <- function(...) {
    deprecate_warn(
        when = "0.2.0",
        what = "spatDimGenePlot3D()",
        with = "spatDimFeatPlot3D()"
    )
    spatDimFeatPlot3D(...)
}
