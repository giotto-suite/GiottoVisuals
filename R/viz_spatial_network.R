#' Bidirectionalize a canonical spatial-network edge table into
#' source/target form and compute `rank_int` (rank-by-distance per source).
#' Replacement for the legacy `convert_to_full_spatial_network()`; the
#' rank is recomputed JIT because any topology mutation (subset, filter)
#' would invalidate a stored value.
#' @keywords internal
#' @noRd
.full_with_rank <- function(sn_dt) {
    distance <- rank_int <- source <- NULL # data.table NSE
    rev <- data.table::copy(sn_dt)
    data.table::setnames(rev, c("from", "to"), c("to", "from"))
    full <- unique(rbind(sn_dt, rev))
    data.table::setnames(full, c("from", "to"), c("source", "target"))
    data.table::setorder(full, source, distance)
    full[, rank_int := seq_len(.N), by = "source"]
    full
}


#' @title Plot spatial distance distribution
#' @name spatNetwDistributionsDistance
#' @description This function return histograms displaying the distance
#' distribution for each spatial k-neighbor
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spatial_network_name name of spatial network
#' @param hist_bins number of binds to use for the histogram
#' @param test_distance_limit effect of different distance threshold on
#' k-neighbors
#' @param ncol number of columns to visualize the histograms in
#' @return ggplot plot
#' @export
spatNetwDistributionsDistance <- function(gobject,
    spat_unit = NULL,
    spatial_network_name = "spatial_network",
    hist_bins = 30,
    test_distance_limit = NULL,
    ncol = 1,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "spatNetwDistributionsDistance") {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )

    # data.table variables
    distance <- rank_int <- status <- label <- keep <- NULL

    ## spatial network
    spatial_network <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "networkDT"
    )

    ## bidirectionalize edges and compute rank_int per source
    spatial_network <- .full_with_rank(spatial_network)

    if (is.null(spatial_network)) {
        stop("spatial network ", spatial_network_name, " was not found")
    }

    if (!is.null(test_distance_limit)) {
        removed_neighbors <- spatial_network[distance > test_distance_limit,
            .N,
            by = rank_int
        ]
        removed_neighbors[, "status" := "remove"]
        keep_neighbors <- spatial_network[distance <= test_distance_limit,
            .N,
            by = rank_int
        ]
        keep_neighbors[, "status" := "keep"]

        dist_removal_dt <- rbind(removed_neighbors, keep_neighbors)
        data.table::setorder(dist_removal_dt, rank_int)

        dist_removal_dt_dcast <- data.table::dcast.data.table(
            data = dist_removal_dt, rank_int ~ status, value.var = "N",
            fill = 0
        )
        dist_removal_dt_dcast[, label := paste0(
            "keep:",
            keep, "\n remove:", remove
        )]
    }

    # text location coordinates
    middle_distance <- max(spatial_network$distance) / (3 / 2)
    freq_dt <- spatial_network[, table(cut(distance, breaks = 30)),
        by = rank_int
    ]
    middle_height <- max(freq_dt$V1) / (3 / 2)

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::labs(title = "distance distribution per k-neighbor")
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::geom_histogram(
        data = spatial_network,
        ggplot2::aes(x = distance),
        color = "white", fill = "black",
        bins = hist_bins
    )
    pl <- pl + ggplot2::facet_wrap(~rank_int, ncol = ncol)
    if (!is.null(test_distance_limit)) {
        pl <- pl + ggplot2::geom_vline(
            xintercept = test_distance_limit,
            color = "red"
        )
        pl <- pl + ggplot2::geom_text(
            data = dist_removal_dt_dcast,
            ggplot2::aes(
                x = middle_distance,
                y = middle_height,
                label = label
            )
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




#' @title spatNetwDistributionsKneighbors
#' @description This function returns a histogram displaying the number
#' of k-neighbors distribution for each cell
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spatial_network_name name of spatial network
#' @param hist_bins number of binds to use for the histogram
#' @return ggplot plot
#' @export
spatNetwDistributionsKneighbors <- function(gobject,
    spat_unit = NULL,
    spatial_network_name = "spatial_network",
    hist_bins = 30,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "spatNetwDistributionsKneighbors") {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )

    # data.table variables
    N <- NULL

    ## spatial network
    # spatial_network = gobject@spatial_network[[spatial_network_name]]
    spatial_network <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "networkDT"
    )

    ## bidirectionalize edges and compute rank_int per source
    spatial_network <- .full_with_rank(spatial_network)

    if (is.null(spatial_network)) {
        stop("spatial network ", spatial_network_name, " was not found")
    }

    spatial_network_dt <- data.table::as.data.table(
        spatial_network[, table(source)]
    )

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::labs(
        title = "k-neighbor distribution for all cells",
        x = "k-neighbors/cell"
    )
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::geom_histogram(
        data = spatial_network_dt,
        ggplot2::aes(x = N), color = "white",
        fill = "black", bins = hist_bins
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



#' @title spatNetwDistributions
#' @name spatNetwDistributions
#' @description This function return histograms displaying the distance
#' distribution for each spatial k-neighbor
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spatial_network_name name of spatial network
#' @param distribution show the distribution of cell-to-cell distance or
#' number of k neighbors
#' @param hist_bins number of binds to use for the histogram
#' @param test_distance_limit effect of different distance threshold on
#' k-neighbors
#' @param ncol number of columns to visualize the histograms in
#' @details The \strong{distance} option shows the spatial distance
#' distribution for each nearest neighbor rank (1st, 2nd, 3th, ... neighbor).
#' With this option the user can also test the effect of a distance limit
#' on the spatial network. This distance limit can be used to remove neighbor
#' cells that are considered to far away. \cr
#' The \strong{k_neighbors} option shows the number of k neighbors distribution
#' over all cells.
#' @return ggplot plot
#' @export
spatNetwDistributions <- function(gobject,
    spat_unit = NULL,
    spatial_network_name = "spatial_network",
    distribution = c("distance", "k_neighbors"),
    hist_bins = 30,
    test_distance_limit = NULL,
    ncol = 1,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "spatNetwDistributions") {
    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )

    ## histogram to show
    distribution <- match.arg(distribution,
        choices = c("distance", "k_neighbors")
    )

    ## spatial network
    spatial_network <- getSpatialNetwork(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spatial_network_name,
        output = "networkDT"
    )
    if (is.null(spatial_network)) {
        stop("spatial network ", spatial_network_name, " was not found")
    }

    switch(distribution,
        "distance" = {
            spatNetwDistributionsDistance(
                gobject = gobject,
                spat_unit = spat_unit,
                spatial_network_name = spatial_network_name,
                hist_bins = hist_bins,
                test_distance_limit = test_distance_limit,
                ncol = ncol,
                show_plot = show_plot,
                return_plot = return_plot,
                save_plot = save_plot,
                save_param = save_param,
                default_save_name = default_save_name
            )
        },
        "k_neighbors" = {
            spatNetwDistributionsKneighbors(
                gobject = gobject,
                spat_unit = spat_unit,
                spatial_network_name = spatial_network_name,
                hist_bins = hist_bins,
                show_plot = show_plot,
                return_plot = return_plot,
                save_plot = save_plot,
                save_param = save_param,
                default_save_name = default_save_name
            )
        }
    )
}






#' @title plotStatDelaunayNetwork
#' @name plotStatDelaunayNetwork
#' @description Plots network statistics for a Delaunay network..
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param method package to use to create a Delaunay network. See
#' `GiottoClass::createSpatialDelaunayNetwork()`.
#' @param dimensions which spatial dimensions to use (maximum 2 dimensions)
#' @param maximum_distance distance cuttof for Delaunay neighbors to consider
#' @param minimum_k minimum neigbhours if maximum_distance != NULL
#' @param options (geometry) String containing extra control options for
#' the underlying Qhull command; see the Qhull documentation
#' (../doc/qhull/html/qdelaun.html) for the available options.
#' (default = 'Pp', do not report precision problems)
#' @param Y (RTriangle) If TRUE prohibits the insertion of Steiner
#' points on the mesh boundary.
#' @param j (RTriangle) If TRUE jettisons vertices that are not part
#' of the final triangulation from the output.
#' @param S (RTriangle) Specifies the maximum number of added Steiner points.
#' @param \dots Other parameters
#' @return giotto object with updated spatial network slot
#' @export
plotStatDelaunayNetwork <- function(gobject,
    feat_type = NULL,
    spat_unit = NULL,
    method = c("delaunayn_geometry", "deldir", "RTriangle"),
    dimensions = "all",
    maximum_distance = "auto", # all
    minimum_k = 0, # all
    options = "Pp", # geometry
    Y = TRUE, # RTriange
    j = TRUE, # RTriange
    S = 0, # RTriange
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "plotStatDelaunayNetwork",
    ...) {
    # data.table variables
    distance <- rank_int <- N <- NULL

    delaunay_network_DT <- createSpatialDelaunayNetwork(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        method = method,
        dimensions = dimensions,
        name = "temp_network",
        maximum_distance = maximum_distance, # all
        minimum_k = minimum_k, # all
        options = options, # geometry
        Y = Y, # RTriange
        j = j, # RTriange
        S = S, # RTriange
        return_gobject = FALSE,
        output = "data.table",
        ...
    )

    delaunay_network_DT_c <- .full_with_rank(delaunay_network_DT)

    ## create visuals
    pl1 <- ggplot2::ggplot(
        delaunay_network_DT,
        ggplot2::aes(x = factor(""), y = distance)
    )
    pl1 <- pl1 + ggplot2::theme_classic() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    pl1 <- pl1 +
        ggplot2::geom_boxplot(outlier.colour = "red", outlier.shape = 1)
    pl1 <- pl1 +
        ggplot2::labs(
            title = "Delaunay network",
            y = "cell-cell distances", x = ""
        )

    pl2 <- ggplot2::ggplot(
        delaunay_network_DT_c,
        ggplot2::aes(x = factor(rank_int), y = distance)
    )
    pl2 <- pl2 + ggplot2::theme_classic() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    pl2 <- pl2 +
        ggplot2::geom_boxplot(outlier.colour = "red", outlier.shape = 1)
    pl2 <- pl2 + ggplot2::labs(
        title = "Delaunay network by neigbor ranking",
        y = "cell-cell distances", x = ""
    )

    neighbors <- delaunay_network_DT_c[, .N, by = source]
    pl3 <- ggplot2::ggplot()
    pl3 <- pl3 + ggplot2::theme_classic() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    pl3 <- pl3 +
        ggplot2::geom_histogram(
            data = neighbors,
            ggplot2::aes(x = as.factor(N)), stat = "count"
        )
    pl3 <- pl3 +
        ggplot2::labs(
            title = "Delaunay network neigbors per cell",
            y = "count", x = ""
        )
    pl3

    savelist <- list(pl1, pl2, pl3)

    ## combine plots with cowplot
    combo_plot <- cowplot::plot_grid(pl1, pl2, NULL, pl3,
        ncol = 2,
        rel_heights = c(1, 1), rel_widths = c(1, 2), align = "v"
    )

    return(plot_output_handler(
        gobject = gobject,
        plot_object = combo_plot,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}
