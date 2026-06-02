## * ####
## 2-D ggplots ####
## ----------- ##



## ** spatial plotting ####





#' @title .spatPlot2D_single
#' @name .spatPlot2D_single
#' @description Visualize cells according to spatial coordinates
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_spatnet_params
#' @inheritParams plot_spatenr_params
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of giotto image(s) to plot
#' @param spat_loc_name name of spatial locations
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparency of voronoi 'cells'
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @param verbose be verbose
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @seealso \code{\link{spatPlot3D}}
.spatPlot2D_single <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        plot_method = "ggplot",
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        spat_loc_name = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        spat_enr_names = NULL,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = "divergent",
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        point_shape = c("border", "no_border", "voronoi"),
        point_size = 3,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        show_cluster_center = FALSE,
        show_center_label = FALSE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        show_network = FALSE,
        spatial_network_name = "Delaunay_network",
        network_color = NULL,
        network_alpha = 1,
        show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        grid_color = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        other_cells_alpha = 0.1,
        coord_fix_ratio = 1,
        title = NULL,
        show_legend = TRUE,
        show_axes = NULL,
        legend_text = 8,
        legend_symbol_size = 1,
        background_color = "white",
        vor_border_color = "white",
        vor_max_radius = 200,
        vor_alpha = 1,
        axis_text = 8,
        axis_title = 8,
        theme_param = list(),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        verbose = FALSE,
        save_param = list(),
        default_save_name = "spatPlot2D_single") {
    # Check params
    .gg_assert_giotto_single(gobject)

    point_shape <- match.arg(
        point_shape,
        choices = c("border", "no_border", "voronoi")
    )

    show_axes <- show_axes %null% TRUE

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

    ## get giotto image(s) ##
    if (isTRUE(show_image) && is.null(gimage)) {
        gimage <- getGiottoImage(
            gobject = gobject,
            name = image_name
        )
    }


    ## get spatial cell locations
    cell_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table",
        copy_obj = TRUE,
        verbose = verbose
    )
    if (is.null(cell_locations)) {
        return(NULL)
    }


    ## extract spatial network
    if (show_network == TRUE) {
        spatial_network <- getSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            name = spatial_network_name,
            output = "networkDT",
            verbose = verbose
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
            name = spatial_grid_name,
            return_grid_Obj = FALSE
        )
    } else {
        spatial_grid <- NULL
    }


    ## get cell metadata

    if (is.null(spat_loc_name)) {
        if (!is.null(slot(gobject, "spatial_locs"))) {
            spat_loc_name <- list_spatial_locations_names(
                gobject,
                spat_unit = spat_unit
            )[[1]]
        } else {
            spat_loc_name <- NULL
            message("No spatial locations have been found")
            return(NULL)
        }
    }

    cell_metadata <- try(
        expr = combineMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            spat_loc_name = spat_loc_name,
            spat_enr_names = spat_enr_names,
            verbose = verbose
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
        cell_locations_metadata_other <-
            cell_locations_metadata[!cell_locations_metadata$cell_ID %in%
                select_cells]
        cell_locations_metadata_selected <-
            cell_locations_metadata[cell_locations_metadata$cell_ID %in%
                select_cells]
        spatial_network <- spatial_network[spatial_network$to %in%
            select_cells & spatial_network$from %in%
            select_cells]

        # if specific cells are selected
        # cell_locations_metadata = cell_locations_metadata_selected
    } else if (is.null(select_cells)) {
        cell_locations_metadata_selected <- cell_locations_metadata
        cell_locations_metadata_other <- NULL
    }


    # update cell_color_code
    # only keep names from selected groups
    if (!is.null(select_cell_groups) & !is.null(cell_color_code)) {
        cell_color_code <- cell_color_code[names(cell_color_code) %in%
            select_cell_groups]
    }

    # data.table and ggplot variables
    sdimx_begin <- sdimy_begin <- sdimx_end <- sdimy_end <- x_start <-
        x_end <- y_start <- y_end <- NULL


    ### create 2D plot with ggplot ###

    if (isTRUE(verbose)) {
        message("Data table with selected information (e.g. cells):")
        message(cell_locations_metadata_selected[seq_len(5), ])

        message("Data table with non-selected information (e.g. cells):")
        message(cell_locations_metadata_other[seq_len(5), ])
    }


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()

    ## plot image ##
    if (isTRUE(show_image) && !is.null(gimage)) {
        pl <- plot_spat_image_layer_ggplot(
            gg_obj = pl,
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            spat_loc_name = spat_loc_name,
            gimage = gimage
        )
    }


    ## plot spatial network
    if (!is.null(spatial_network) && isTRUE(show_network)) {
        if (is.null(network_color)) network_color <- "red"
        pl <- pl + ggplot2::geom_segment(
            data = spatial_network,
            aes(
                x = sdimx_begin,
                y = sdimy_begin,
                xend = sdimx_end,
                yend = sdimy_end
            ),
            color = network_color,
            linewidth = 0.5,
            alpha = network_alpha
        )
    }


    ## plot spatial grid
    if (!is.null(spatial_grid) && isTRUE(show_grid)) {
        if (is.null(grid_color)) grid_color <- "black"
        pl <- pl + ggplot2::geom_rect(
            data = spatial_grid,
            aes(
                xmin = x_start,
                xmax = x_end,
                ymin = y_start,
                ymax = y_end
            ),
            color = grid_color,
            fill = NA
        )
    }


    ## plot point layer
    point_general_params <- list(
        ggobject = pl,
        ext = ext(gobject, prefer = "spatlocs"),
        instrs = instructions(gobject),
        sdimx = sdimx,
        sdimy = sdimy,
        plot_method = plot_method,
        cell_locations_metadata_selected = cell_locations_metadata_selected,
        cell_locations_metadata_other = cell_locations_metadata_other,
        cell_color = cell_color,
        color_as_factor = color_as_factor,
        cell_color_code = cell_color_code,
        cell_color_gradient = cell_color_gradient,
        gradient_midpoint = gradient_midpoint,
        gradient_style = gradient_style,
        gradient_limits = gradient_limits,
        select_cell_groups = select_cell_groups,
        select_cells = select_cells,
        point_size = point_size,
        point_alpha = point_alpha,
        show_cluster_center = show_cluster_center,
        show_center_label = show_center_label,
        center_point_size = center_point_size,
        label_size = label_size,
        label_fontface = label_fontface,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        other_point_size = other_point_size,
        show_legend = show_legend
    )

    point_border_specific_params <- list(
        point_border_stroke = point_border_stroke, # specific
        point_border_col = point_border_col, # specific
        center_point_border_col = center_point_border_col, # specific
        center_point_border_stroke = center_point_border_stroke # specific
    )

    point_voronoi_specific_params <- list(
        background_color = background_color, # specific
        vor_border_color = vor_border_color, # specific
        vor_max_radius = vor_max_radius, # specific
        vor_alpha = vor_alpha # specific
    )

    pl <- switch(point_shape,
        "border" = do.call(
            plot_spat_point_layer_ggplot,
            args = c(
                point_general_params,
                point_border_specific_params
            )
        ),
        "no_border" = do.call(
            plot_spat_point_layer_ggplot_noFILL,
            args = point_general_params
        ),
        "voronoi" = do.call(
            plot_spat_voronoi_layer_ggplot,
            args = c(
                point_general_params,
                point_voronoi_specific_params
            )
        )
    )


    ## adjust theme settings
    gg_theme_args <- c(
        theme_param,
        legend_text = legend_text,
        axis_title = axis_title,
        axis_text = axis_text,
        background_color = background_color
    )
    pl <- pl + do.call(.gg_theme, args = gg_theme_args)

    if (!show_axes) {
        pl <- .theme_remove_axes(pl)
    }

    ## change symbol size of legend
    if (isTRUE(color_as_factor)) {
        if (point_shape %in% c("border", "voronoi")) {
            pl <- pl +
                guides(fill = guide_legend(
                    override.aes = list(size = legend_symbol_size)
                ))
        } else if (point_shape == "no_border") {
            pl <- pl +
                guides(color = guide_legend(
                    override.aes = list(size = legend_symbol_size)
                ))
        }
    }


    # fix coord ratio
    pl <- .aspect_ratio(pl, coord_fix_ratio)

    # provide x, y and plot titles
    if (is.null(title)) title <- cell_color
    pl <- pl + ggplot2::labs(
        x = "x coordinates", y = "y coordinates",
        title = title
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







#' @rdname spatPlot
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_image_params
#' @inheritParams plot_spatenr_params
#' @inheritParams plot_params
#' @param plot_method method to plot points. Either "ggplot" (default) or
#' "scattermore" (rasterized and faster for large datasets)
#' @param spat_loc_name name of spatial locations
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on
#' cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_shape shape of points (border, no_border or voronoi)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of point
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param center_point_border_col border color of center points
#' @param center_point_border_stroke border stroke size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_network show underlying spatial network
#' @param spatial_network_name name of spatial network to use
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size point size of not selected cells
#' @param other_cells_alpha alpha of not selected cells
#' @param coord_fix_ratio fix ratio between x and y-axis (default = 1)
#' @param title title of plot
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param vor_border_color border color for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparency of voronoi 'cells'
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @param view,space optional [GiottoClass::giottoView-class] /
#' [GiottoClass::giottoSpace-class] or the name of one slotted on
#' `gobject`. When supplied, `gobject` is pre-narrowed via
#' [GiottoClass::materialize()] before data is fetched. One resolver
#' pass shared across all slot reads; same code path for in-mem and
#' on-disk gobjects.
#' @details coord_fix_ratio: set to NULL to use default ggplot parameters
#' @returns ggplot
#' @export
spatPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        plot_method = "ggplot",
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        group_by = NULL,
        group_by_subset = NULL,
        spat_loc_name = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        spat_enr_names = NULL,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        point_shape = c("border", "no_border", "voronoi"),
        point_size = 3,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        show_cluster_center = FALSE,
        show_center_label = FALSE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        show_network = FALSE,
        spatial_network_name = "Delaunay_network",
        network_color = NULL,
        network_alpha = 1,
        show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        grid_color = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        other_cells_alpha = 0.1,
        coord_fix_ratio = 1,
        title = NULL,
        show_legend = TRUE,
        show_axes = NULL,
        legend_text = 10,
        legend_symbol_size = 2,
        background_color = "white",
        vor_border_color = "white",
        vor_max_radius = 200,
        vor_alpha = 1,
        axis_text = 8,
        axis_title = 8,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        theme_param = list(),
        default_save_name = "spatPlot2D",
        view = NULL,
        space = NULL) {
    if (inherits(gobject, "giottoMulti")) {
        return(.gg_multi_dispatch_spatial(
            plot_fn = spatPlot2D,
            named = mget(names(formals())),
            gobject = gobject, view = view, space = space
        ))
    }
    .gg_assert_giotto_single(gobject)

    # Pre-narrow once for the slots this plot reads.
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "dimension_reduction",
            "images"))

    # deprecation message
    if (!is.null(largeImage_name)) {
        deprecate_warn(
            when = "0.2.0",
            what = "spatPlot2D(largeImage_name)",
            details = c(
                "Use `image_name` argument instead for all images to plot."
            )
        )
        image_name <- c(image_name, largeImage_name)
    }

    point_shape <- match.arg(
        point_shape, c("border", "no_border", "voronoi")
    )

    if (identical(plot_method, "scattermore") &&
        point_shape != "no_border") {
        warning("point_shape changed to \"no_border\" for scattermore")
        point_shape <- "no_border"
    }

    # create args list needed for each call to .spatPlot2D_single()
    # 1. - grab all params available
    # 2. - subset to those needed
    spp_params <- get_args_list(keep = c(
        # [gobject params]
        "gobject", "spat_unit", "feat_type",
        # [image params]
        "show_image", "gimage", "image_name",
        # [spatlocs params]
        "spat_loc_name", "sdimx", "sdimy",
        # [access spatial enrichments]
        "spat_enr_names",
        # [point aes]
        "plot_method", "cell_color", "color_as_factor", "cell_color_code",
        "cell_color_gradient",
        "gradient_midpoint", "gradient_style", "gradient_limits",
        "point_shape", "point_size", "point_alpha", "point_border_col",
        "point_border_stroke",
        # [select cell params]
        "select_cell_groups", "select_cells",
        # [voronoi-point params]
        "vor_border_color", "vor_max_radius", "vor_alpha",
        # [others aes]
        "show_other_cells", "other_cell_color", "other_point_size",
        "other_cells_alpha",
        # [cluster aes]
        "show_cluster_center", "show_center_label", "center_point_size",
        "center_point_border_col", "center_point_border_stroke",
        # [label aes]
        "label_size", "label_fontface",
        # [network aes]
        "show_network", "spatial_network_name", "network_color",
        "network_alpha",
        # [grid aes]
        "show_grid", "spatial_grid_name", "grid_color",
        # [figure params]
        "coord_fix_ratio", "show_legend", "show_axes", "legend_text",
        "legend_symbol_size", "background_color", "axis_text",
        "axis_title", "title",
        # [return params]
        "show_plot", "return_plot", "save_plot", "save_param",
        "default_save_name",
        # [gg params]
        "theme_param"
    ))


    ## check group_by
    if (is.null(group_by)) { # ----------------------------------------------- #

        handle_errors(
            do.call(.spatPlot2D_single, args = spp_params),
            prefix = "spatPlot2D()"
        )

    } else { # -------------------------------------------------------------- #

        # setup for group_by
        # params relevant for plotting that are updated in this section prior
        # to the for loop MUST be updated in group_by static settings section

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
        # ! update spat_unit & feat_type in static settings !

        ## check metadata for valid group_by information
        comb_metadata <- combineMetadata(
            gobject = gobject,
            spat_loc_name = spat_loc_name,
            feat_type = feat_type,
            spat_unit = spat_unit,
            spat_enr_names = spat_enr_names
        )
        possible_meta_groups <- colnames(comb_metadata)

        ## error if group_by col is not found
        if (!group_by %in% possible_meta_groups) {
            stop("group_by ", group_by, " was not found in pDataDT()")
        }

        unique_groups <- unique(comb_metadata[[group_by]])

        # subset unique_groups
        # These unique_groups will be used to iterate through subsetting then
        # plotting the giotto object multiple times.
        if (!is.null(group_by_subset)) {
            not_found <- group_by_subset[!group_by_subset %in% unique_groups]
            if (length(not_found) > 0) {
                message("the following subset was not found: ", not_found)
            }
            unique_groups <- unique_groups[unique_groups %in% group_by_subset]
        }

        # create matching cell_color_code
        if (is.null(cell_color_code)) {
            if (is.character(cell_color)) {
                if (cell_color %in% colnames(comb_metadata)) {
                    if (isTRUE(color_as_factor)) {
                        number_colors <- length(
                            unique(comb_metadata[[cell_color]])
                        )
                        cell_color_code <- set_default_color_discrete_cell(
                            instrs = instructions(gobject)
                        )(n = number_colors)
                        names(cell_color_code) <- unique(
                            comb_metadata[[cell_color]]
                        )
                        cell_color_code <- cell_color_code
                    }
                }
            }
        }
        # ! update cell_color_code in static settings !



        ## plotting ##
        savelist <- list()

        # group_by static settings #
        # update these params
        spp_params$spat_unit <- spat_unit
        spp_params$feat_type <- feat_type
        spp_params$cell_color_code <- cell_color_code
        # apply group_by specific defaults
        spp_params$show_plot <- FALSE
        spp_params$return_plot <- TRUE
        spp_params$save_plot <- FALSE
        spp_params$save_param <- list()
        spp_params$default_save_name <- "spatPlot2D"


        for (group_id in seq_along(unique_groups)) {
            group <- unique_groups[group_id]

            subset_cell_IDs <- comb_metadata[
                get(group_by) == group][["cell_ID"]]
            spp_params$gobject <- subsetGiotto(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                cell_ids = subset_cell_IDs,
                verbose = FALSE
            )

            # use a different image per group if there are the same
            # number of names provided as there are groups
            # Otherwise, use the same image (or NULL) for all groups (default)
            if (length(unique_groups) == length(image_name)) {
                spp_params$image_name <- image_name[group_id]
            }

            pl <- handle_errors(
                do.call(.spatPlot2D_single, args = spp_params),
                prefix = sprintf("spatPlot2D() group: %s", group_id)
            )

            savelist[[group_id]] <- pl
        }

        # combine plots with cowplot
        combo_plot <- cowplot::plot_grid(
            plotlist = savelist,
            ncol = set_default_cow_n_col(
                cow_n_col = cow_n_col,
                nr_plots = length(savelist)
            ),
            rel_heights = cow_rel_h,
            rel_widths = cow_rel_w,
            align = cow_align
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
    } # --------------------------------------------------------------------- #
}




#' @title spatPlot
#' @name spatPlot
#' @description Visualize cells according to spatial coordinates
#' @param \dots spatPLot(...) passes to spatPlot2D
#' @return ggplot (2D), plotly (3D)
#' @family spatial visualizations
#' @returns ggplot
#' @examples
#' # 2D Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatPlot(g, show_image = TRUE, image_name = "image")
#'
#' # the more specific spatPlot2D with networks shown
#' spatPlot2D(g, show_image = TRUE, image_name = "image", show_network = TRUE)
#'
#' # plotting of some cell metadata (number of different features detected)
#' spatPlot2D(g,
#'     show_image = TRUE,
#'     image_name = "image",
#'     cell_color = "nr_feats",
#'     color_as_factor = FALSE,
#'     gradient_style = "sequential"
#' )
#'
#'
#' # 3D Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' starmap <- GiottoData::loadGiottoMini("starmap", verbose = FALSE)
#'
#' # default is to rescale plot as a 3D cube
#' spatPlot3D(starmap, cell_color = "leiden_clus")
#' # real scaling
#' spatPlot3D(starmap, cell_color = "leiden_clus", axis_scale = "real")
#'
#' # plot with selected cell groups
#' spatPlot3D(starmap,
#'     cell_color = "cell_types",
#'     color_as_factor = TRUE,
#'     select_cell_groups = c("cell F", "cell C", "cell A"),
#'     other_point_size = 1
#' )
#'
#' # use the "sequential" style gradient default
#' spatPlot3D(starmap,
#'     cell_color = "total_expr",
#'     color_as_factor = FALSE,
#'     point_alpha = 0.5,
#'     axis_scale = "real",
#'     gradient_style = "sequential"
#' )
#'
#' # specific color gradient
#' spatPlot3D(starmap,
#'     cell_color = "total_expr",
#'     color_as_factor = FALSE,
#'     point_alpha = 0.7,
#'     axis_scale = "cube",
#'     cell_color_gradient = "mako"
#' )
#' @export
#' @seealso \code{\link{spatPlot3D}}
spatPlot <- function(...) {
    spatPlot2D(...)
}








## ** spatial deconvolution plotting ####


#' @title spatDeconvPlot
#' @name spatDeconvPlot
#' @description Visualize cell type enrichment / deconvolution results
#' in a scatterpie
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @param deconv_name name of deconvolution results to use
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param spat_loc_name name of spatial locations
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param line_color color of line within pie charts
#' @param pie_scale amount to scale the pie size if there is no radius
#' mapping exists
#' @param alpha alpha of pie charts
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param legend_text size of legend text
#' @param background_color color of plot background
#' @param title title for plot (default = deconv_name)
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @returns ggplot
#' @export
spatDeconvPlot <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        deconv_name = "DWLS",
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        spat_loc_name = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        cell_color_code = NULL,
        line_color = NA,
        pie_scale = 0.35,
        alpha = 1,
        legend_text = 8,
        background_color = "white",
        title = NULL,
        axis_text = 8,
        axis_title = 8,
        coord_fix_ratio = 1,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        theme_param = list(),
        default_save_name = "spatDeconvPlot",
        view = NULL,
        space = NULL) {
    if (inherits(gobject, "giottoMulti")) {
        return(.gg_multi_dispatch_spatial(
            plot_fn = spatDeconvPlot,
            named = mget(names(formals())),
            gobject = gobject, view = view, space = space
        ))
    }
    .gg_assert_giotto_single(gobject)
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "images"))
    # check for installed packages
    package_check(pkg_name = "scatterpie", repository = "CRAN")

    # deprecation message
    if (!is.null(largeImage_name)) {
        deprecate_warn(
            when = "0.2.0",
            what = "spatDeconvPlot(largeImage_name)",
            details = c(
                "Use `image_name` argument instead for all images to plot."
            )
        )
        image_name <- c(image_name, largeImage_name)
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

    ## get giotto image(s) ##
    if (isTRUE(show_image) && is.null(gimage)) {
        gimage <- getGiottoImage(
            gobject = gobject,
            name = image_name
        )
    }


    ## get spatial cell locations
    spatial_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table"
    )
    if (is.null(spatial_locations)) {
        return(NULL)
    }

    ## deconvolution results
    spatial_enrichment <- getSpatialEnrichment(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        name = deconv_name,
        output = "data.table"
    )

    ### create 2D plot with ggplot ###

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()

    ## plot image ##
    if (isTRUE(show_image) && !is.null(gimage)) {
        pl <- plot_spat_image_layer_ggplot(
            gg_obj = pl,
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            spat_loc_name = spat_loc_name,
            gimage = gimage
        )
    }


    ## plot scatterpie ##
    pl <- plot_spat_scatterpie_layer_ggplot(
        ggobject = pl,
        instrs = instructions(gobject),
        sdimx = sdimx,
        sdimy = sdimy,
        spatial_locations = spatial_locations,
        spatial_enrichment = spatial_enrichment,
        pie_scale = pie_scale,
        color = line_color,
        alpha = alpha,
        cell_color_code = cell_color_code
    )


    ## adjust theme setting
    gg_theme_args <- c(
        theme_param,
        legend_text = legend_text,
        axis_title = axis_title,
        axis_text = axis_text,
        background_color = background_color
    )
    pl <- pl + do.call(.gg_theme, args = gg_theme_args)

    # fix coord ratio
    pl <- .aspect_ratio(pl, coord_fix_ratio)

    # provide x, y and plot titles
    if (is.null(title)) title <- deconv_name
    pl <- pl +
        ggplot2::labs(x = "x coordinates", y = "y coordinates", title = title)


    # print, return and save parameters
    show_plot <- ifelse(is.null(show_plot),
        instructions(gobject, param = "show_plot"),
        show_plot
    )
    save_plot <- ifelse(is.null(save_plot),
        instructions(gobject, param = "save_plot"),
        save_plot
    )
    return_plot <- ifelse(is.null(return_plot),
        instructions(gobject, param = "return_plot"),
        return_plot
    )

    ## print plot
    if (show_plot == TRUE) {
        print(pl)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call("all_plots_save_function", c(list(
            gobject = gobject,
            plot_object = pl,
            default_save_name = default_save_name
        ), save_param))
    }

    ## return plot
    if (return_plot == TRUE) {
        return(pl)
    }
}





# ** dim reduction plotting ####




# Create a single 2D dimplot. This is looped through by dimPlot2D() if needed.
#' @noRd
#' @keywords internal
.dimPlot2D_single <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = NULL,
        dim1_to_use = 1,
        dim2_to_use = 2,
        spat_enr_names = NULL,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 0.5,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        edge_alpha = NULL,
        point_shape = c("border", "no_border"),
        point_size = 1,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        title = NULL,
        show_legend = TRUE,
        legend_text = 8,
        legend_symbol_size = 1,
        background_color = "white",
        axis_text = 8,
        axis_title = 8,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "dimPlot2D_single") {
    # giottoMulti pass-through: dim reductions live on the joint
    # @dimension_reduction slot (built by runPCA/runUMAP on the
    # assembled expression). All downstream getters return joint
    # subobjects — ggplot sees one combined embedding.

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

    # specify dim_reduction_name according to provided feat_type
    if (!is.null(dim_reduction_to_use)) {
        if (is.null(dim_reduction_name)) {
            if (feat_type == "rna") {
                dim_reduction_name <- dim_reduction_to_use
            } else {
                dim_reduction_name <- paste0(
                    feat_type, ".",
                    dim_reduction_to_use
                )
            }
        }
    }

    ## point shape ##
    point_shape <- match.arg(point_shape, c("border", "no_border"))

    ## dimension reduction ##
    # test if dimension reduction was performed

    dim_red_names <- list_dim_reductions_names(
        gobject = gobject, data_type = "cells",
        spat_unit = spat_unit, feat_type = feat_type,
        dim_type = dim_reduction_to_use
    )

    if (!dim_reduction_name %in% dim_red_names) {
        stop(
            "\n dimension reduction: ", dim_reduction_to_use,
            " or dimension reduction name: ", dim_reduction_name,
            " is not available \n"
        )
    }


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

    # data.table variables
    cell_ID <- NULL

    dim_DT <- data.table::as.data.table(dim_dfr)
    dim_DT[, cell_ID := as.character(rownames(dim_dfr))]

    ## annotated cell metadata
    cell_metadata <- combineMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spat_enr_names = spat_enr_names,
        spat_loc_name = NULL
    )

    cell_metadata[, cell_ID := as.character(cell_ID)]

    annotated_DT <- data.table::merge.data.table(cell_metadata,
        dim_DT,
        by = "cell_ID"
    )


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

        network_DT <- data.table::as.data.table(
            igraph::as_data_frame(selected_nn_network, what = "edges")
        )

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

    # add % variance information if reduction is PCA
    if (dim_reduction_to_use == "pca") {
        pcaObj <- getDimReduction(gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "dimObj"
        )
        eigenvalues <- pcaObj@misc$eigenvalues

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
        annotated_DT <- annotated_DT_selected
    }

    ## if no subsets are required
    if (is.null(select_cells) & is.null(select_cell_groups)) {
        annotated_DT_selected <- annotated_DT
        annotated_DT_other <- NULL
    }



    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    ## add network layer
    if (show_NN_network == TRUE) {
        pl <- plot_network_layer_ggplot(
            ggobject = pl,
            instrs = instructions(gobject),
            annotated_network_DT = annotated_network_DT,
            edge_alpha = edge_alpha,
            show_legend = show_legend
        )
    }

    # return(list(pl, annotated_DT_selected, annotated_DT_other))

    if (point_shape == "border") {
        ## add point layer
        pl <- plot_point_layer_ggplot(
            ggobject = pl,
            instrs = instructions(gobject),
            annotated_DT_selected = annotated_DT_selected,
            annotated_DT_other = annotated_DT_other,
            cell_color = cell_color,
            color_as_factor = color_as_factor,
            cell_color_code = cell_color_code,
            cell_color_gradient = cell_color_gradient,
            gradient_midpoint = gradient_midpoint,
            gradient_style = gradient_style,
            gradient_limits = gradient_limits,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_cluster_center = show_cluster_center,
            show_center_label = show_center_label,
            center_point_size = center_point_size,
            center_point_border_col = center_point_border_col,
            center_point_border_stroke = center_point_border_stroke,
            label_size = label_size,
            label_fontface = label_fontface,
            edge_alpha = edge_alpha,
            point_size = point_size,
            point_alpha = point_alpha,
            point_border_col = point_border_col,
            point_border_stroke = point_border_stroke,
            show_legend = show_legend
        )
    } else if (point_shape == "no_border") {
        pl <- plot_point_layer_ggplot_noFILL(
            ggobject = pl,
            instrs = instructions(gobject),
            annotated_DT_selected = annotated_DT_selected,
            annotated_DT_other = annotated_DT_other,
            cell_color = cell_color,
            color_as_factor = color_as_factor,
            cell_color_code = cell_color_code,
            cell_color_gradient = cell_color_gradient,
            gradient_midpoint = gradient_midpoint,
            gradient_style = gradient_style,
            gradient_limits = gradient_limits,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_cluster_center = show_cluster_center,
            show_center_label = show_center_label,
            center_point_size = center_point_size,
            label_size = label_size,
            label_fontface = label_fontface,
            edge_alpha = edge_alpha,
            point_size = point_size,
            point_alpha = point_alpha,
            show_legend = show_legend
        )
    }


    ## add % variance explained to names of plot for PCA ##
    if (dim_reduction_to_use == "pca") {
        if (!is.null(eigenvalues)) {
            x_name <- paste0("pca", "-", dim_names[1])
            y_name <- paste0("pca", "-", dim_names[2])

            # provide x, y and plot titles
            x_title <- sprintf(
                "%s explains %.02f%% of variance",
                x_name, var_expl_vec[dim1_to_use]
            )
            y_title <- sprintf(
                "%s explains %.02f%% of variance",
                y_name, var_expl_vec[dim2_to_use]
            )

            if (is.null(title)) title <- cell_color
            pl <- pl + ggplot2::labs(x = x_title, y = y_title, title = title)
        }
    } else {
        # provide x, y and plot titles
        x_title <- paste0(dim_reduction_to_use, "-", dim_names[1])
        y_title <- paste0(dim_reduction_to_use, "-", dim_names[2])

        if (is.null(title)) title <- cell_color
        pl <- pl + ggplot2::labs(x = x_title, y = y_title, title = title)
    }

    ## adjust titles
    pl <- pl + ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_text),
        axis.text = element_text(size = axis_text),
        axis.title = element_text(size = axis_title),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = background_color)
    )

    ## change symbol size of legend
    if (color_as_factor == TRUE) {
        if (point_shape == "border") {
            pl <- pl + guides(fill = guide_legend(
                override.aes = list(size = legend_symbol_size)
            ))
        } else if (point_shape == "no_border") {
            pl <- pl + guides(color = guide_legend(
                override.aes = list(size = legend_symbol_size)
            ))
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




#' @rdname dimPlot
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_dimred_params
#' @inheritParams plot_nn_net_params
#' @inheritParams plot_spatenr_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_params
#' @returns ggplot
#' @family reduced dimension visualizations
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' dimPlot2D(g)
#' @export
dimPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        group_by = NULL,
        group_by_subset = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = NULL,
        dim1_to_use = 1,
        dim2_to_use = 2,
        spat_enr_names = NULL,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 0.5,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        edge_alpha = NULL,
        point_shape = c("border", "no_border"),
        point_size = 1,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        title = NULL,
        show_legend = TRUE,
        legend_text = 10,
        legend_symbol_size = 2,
        background_color = "white",
        axis_text = 8,
        axis_title = 8,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "dimPlot2D",
        view = NULL,
        space = NULL) {
    # giottoMulti pass-through: gmulti carries a joint dim_reduction
    # slot (PCA / UMAP / TSNE run on the assembled expression matrix);
    # getDimReduction(mg, ...) returns that single joint embedding and
    # downstream getters return joint subobjects. No per-sample loop.
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "dimension_reduction",
            "spatial_enrichment", "expression"))

    handle_errors({
        ## check group_by
        if (is.null(group_by)) {
            .dimPlot2D_single(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                dim_reduction_to_use = dim_reduction_to_use,
                dim_reduction_name = dim_reduction_name,
                dim1_to_use = dim1_to_use,
                dim2_to_use = dim2_to_use,
                spat_enr_names = spat_enr_names,
                show_NN_network = show_NN_network,
                nn_network_to_use = nn_network_to_use,
                network_name = network_name,
                cell_color = cell_color,
                color_as_factor = color_as_factor,
                cell_color_code = cell_color_code,
                cell_color_gradient = cell_color_gradient,
                gradient_midpoint = gradient_midpoint,
                gradient_style = gradient_style,
                gradient_limits = gradient_limits,
                select_cell_groups = select_cell_groups,
                select_cells = select_cells,
                show_other_cells = show_other_cells,
                other_cell_color = other_cell_color,
                other_point_size = other_point_size,
                show_cluster_center = show_cluster_center,
                show_center_label = show_center_label,
                center_point_size = center_point_size,
                center_point_border_col = center_point_border_col,
                center_point_border_stroke = center_point_border_stroke,
                label_size = label_size,
                label_fontface = label_fontface,
                edge_alpha = edge_alpha,
                point_shape = point_shape,
                point_size = point_size,
                point_alpha = point_alpha,
                point_border_col = point_border_col,
                point_border_stroke = point_border_stroke,
                title = title,
                show_legend = show_legend,
                legend_text = legend_text,
                legend_symbol_size = legend_symbol_size,
                background_color = background_color,
                axis_text = axis_text,
                axis_title = axis_title,
                show_plot = show_plot,
                return_plot = return_plot,
                save_plot = save_plot,
                save_param = save_param,
                default_save_name = default_save_name
            )
        } else {
            comb_metadata <- combineMetadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                spat_enr_names = spat_enr_names,
                spat_loc_name = NULL
            )
            possible_meta_groups <- colnames(comb_metadata)

            ## check if group_by is found
            if (!group_by %in% possible_meta_groups) {
                stop("group_by ", group_by, " was not found in pDataDT()")
            }

            unique_groups <- unique(comb_metadata[[group_by]])

            # subset unique_groups
            if (!is.null(group_by_subset)) {
                not_found <- group_by_subset[
                    !group_by_subset %in% unique_groups]

                if (length(not_found) > 0) {
                    message("the following subset was not found: ", not_found)
                }
                unique_groups <- unique_groups[
                    unique_groups %in% group_by_subset]
            }


            # create matching cell_color_code for groupby factors
            # best done prior to the following groupby subsetGiotto() operation
            if (is.null(cell_color_code)) { # TODO add getColors() support
                if (is.character(cell_color)) {
                    if (cell_color %in% colnames(comb_metadata)) {
                        if (color_as_factor == TRUE) {
                            number_colors <- length(
                                unique(comb_metadata[[cell_color]])
                            )
                            cell_color_code <- set_default_color_discrete_cell(
                                instrs = instructions(gobject)
                            )(n = number_colors)
                            names(cell_color_code) <- unique(
                                comb_metadata[[cell_color]]
                            )
                            cell_color_code <- cell_color_code
                        }
                    }
                }
            }

            ## plotting ##
            savelist <- list()


            for (group_id in seq_len(length(unique_groups))) {
                group <- unique_groups[group_id]

                subset_cell_IDs <- comb_metadata[
                    get(group_by) == group
                ][["cell_ID"]]
                temp_gobject <- subsetGiotto(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    cell_ids = subset_cell_IDs
                )

                pl <- .dimPlot2D_single(
                    gobject = temp_gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    dim_reduction_to_use = dim_reduction_to_use,
                    dim_reduction_name = dim_reduction_name,
                    dim1_to_use = dim1_to_use,
                    dim2_to_use = dim2_to_use,
                    spat_enr_names = spat_enr_names,
                    show_NN_network = show_NN_network,
                    nn_network_to_use = nn_network_to_use,
                    network_name = network_name,
                    cell_color = cell_color,
                    cell_color_code = cell_color_code,
                    color_as_factor = color_as_factor,
                    cell_color_gradient = cell_color_gradient,
                    gradient_midpoint = gradient_midpoint,
                    gradient_style = gradient_style,
                    gradient_limits = gradient_limits,
                    select_cell_groups = select_cell_groups,
                    select_cells = select_cells,
                    show_other_cells = show_other_cells,
                    other_cell_color = other_cell_color,
                    other_point_size = other_point_size,
                    show_cluster_center = show_cluster_center,
                    show_center_label = show_center_label,
                    center_point_size = center_point_size,
                    center_point_border_col = center_point_border_col,
                    center_point_border_stroke = center_point_border_stroke,
                    label_size = label_size,
                    label_fontface = label_fontface,
                    edge_alpha = edge_alpha,
                    point_shape = point_shape,
                    point_size = point_size,
                    point_alpha = point_alpha,
                    point_border_col = point_border_col,
                    point_border_stroke = point_border_stroke,
                    title = group,
                    show_legend = show_legend,
                    legend_text = legend_text,
                    legend_symbol_size = legend_symbol_size,
                    background_color = background_color,
                    axis_text = axis_text,
                    axis_title = axis_title,
                    show_plot = FALSE,
                    return_plot = TRUE,
                    save_plot = FALSE,
                    save_param = list(),
                    default_save_name = default_save_name
                )


                savelist[[group_id]] <- pl
            }

            # combine plots with cowplot
            combo_plot <- cowplot::plot_grid(
                plotlist = savelist,
                ncol = set_default_cow_n_col(
                    cow_n_col = cow_n_col,
                    nr_plots = length(savelist)
                ),
                rel_heights = cow_rel_h,
                rel_widths = cow_rel_w,
                align = cow_align
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
    }, prefix = "dimPlot2D()")
}





#' @title Plot dimension reduction
#' @name dimPlot
#' @param \dots dimPlot(...) passes to dimPlot2D()
#' @description Visualize cells according to dimension reduction coordinates
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' dimPlot(g)
#' @export
dimPlot <- function(...) {
    dimPlot2D(...)
}






#' @title plotUMAP_2D
#' @name plotUMAP_2D
#' @description Short wrapper for UMAP visualization
#' @inheritParams data_access_params
#' @param dim_reduction_name name of UMAP
#' @param default_save_name default save name of UMAP plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters, see \code{\link{dimPlot2D}}.
#' For 3D plots see \code{\link{plotUMAP_3D}}
#' @family reduced dimension visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotUMAP_2D(g)
#' @export
plotUMAP_2D <- function(
        gobject,
        dim_reduction_name = NULL,
        default_save_name = "UMAP_2D",
        ...) {
    dimPlot2D(
        gobject = gobject,
        dim_reduction_to_use = "umap",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}


#' @title plotUMAP
#' @name plotUMAP
#' @description Short wrapper for UMAP visualization
#' @inheritParams data_access_params
#' @param dim_reduction_name name of UMAP
#' @param default_save_name default save name of UMAP plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @family reduced dimension visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotUMAP(g)
#'
#' @export
plotUMAP <- function(
        gobject,
        dim_reduction_name = NULL,
        default_save_name = "UMAP",
        ...) {
    dimPlot2D(
        gobject = gobject,
        dim_reduction_to_use = "umap",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}





#' @title plotTSNE_2D
#' @name plotTSNE_2D
#' @description Short wrapper for tSNE visualization
#' @inheritParams data_access_params
#' @param dim_reduction_name name of TSNE
#' @param default_save_name default save name of TSNE plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters, see \code{\link{dimPlot2D}}.
#' For 3D plots see \code{\link{plotTSNE_3D}}
#' @family reduced dimension visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotTSNE_2D(g)
#'
#' @export
plotTSNE_2D <- function(
        gobject,
        dim_reduction_name = NULL,
        default_save_name = "tSNE_2D",
        ...) {
    dimPlot2D(
        gobject = gobject,
        dim_reduction_to_use = "tsne",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}

#' @title plotTSNE
#' @name plotTSNE
#' @description Short wrapper for tSNE visualization
#' @inheritParams data_access_params
#' @param dim_reduction_name name of TSNE
#' @param default_save_name default save name of TSNE plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters, see \code{\link{dimPlot2D}}.
#' For 3D plots see \code{\link{plotTSNE_3D}}
#' @family reduced dimension visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotTSNE(g)
#'
#' @export
plotTSNE <- function(
        gobject,
        dim_reduction_name = NULL,
        default_save_name = "tSNE",
        ...) {
    dimPlot2D(
        gobject = gobject,
        dim_reduction_to_use = "tsne",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}



#' @title plotPCA_2D
#' @name plotPCA_2D
#' @description Short wrapper for PCA visualization
#' @inheritParams data_access_params
#' @param dim_reduction_name name of PCA
#' @param default_save_name default save name of PCA plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters, see \code{\link{dimPlot2D}}.
#' For 3D plots see \code{\link{plotPCA_3D}}
#' @family reduced dimension visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotPCA_2D(g)
#'
#' @export
plotPCA_2D <- function(
        gobject,
        dim_reduction_name = NULL,
        default_save_name = "PCA_2D",
        ...) {
    dimPlot2D(
        gobject = gobject,
        dim_reduction_to_use = "pca",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}



#' @title plotPCA
#' @name plotPCA
#' @description Short wrapper for PCA visualization
#' @inheritParams data_access_params
#' @param dim_reduction_name name of PCA
#' @param default_save_name default save name of PCA plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use
#' -dim_reduction_name -default_save_name
#' @details Description of parameters, see \code{\link{dimPlot2D}}.
#' For 3D plots see \code{\link{plotPCA_3D}}
#' @family reduced dimension visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' plotPCA(g)
#'
#' @export
plotPCA <- function(
        gobject,
        dim_reduction_name = NULL,
        default_save_name = "PCA",
        ...) {
    dimPlot2D(
        gobject = gobject,
        dim_reduction_to_use = "pca",
        dim_reduction_name = dim_reduction_name,
        default_save_name = default_save_name,
        ...
    )
}












## ** spatial and dim reduction plotting ####


#' @title spatDimPlot
#' @name spatDimPlot
#' @description Visualize cells according to spatial AND dimension reduction
#' coordinates 2D
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_spatenr_params
#' @inheritParams plot_dimred_params
#' @inheritParams plot_nn_net_params
#' @inheritParams plot_spatnet_params
#' @inheritParams plot_image_params
#' @inheritParams plot_params
#' @param largeImage_name deprecated
#' @param spat_loc_name name of spatial locations
#' @param plot_alignment direction to align plot
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param spat_point_shape shape of points (border, no_border or voronoi)
#' @param spat_point_size size of spatial points
#' @param spat_point_alpha transparency of spatial points
#' @param spat_point_border_col border color of spatial points
#' @param spat_point_border_stroke border stroke of spatial points
#' @param dim_show_cluster_center show the center of each cluster
#' @param dim_show_center_label provide a label for each cluster
#' @param dim_center_point_size size of the center point
#' @param dim_center_point_border_col border color of center point
#' @param dim_center_point_border_stroke stroke size of center point
#' @param dim_label_size size of the center label
#' @param dim_label_fontface font of the center label
#' @param spat_show_cluster_center show the center of each cluster
#' @param spat_show_center_label provide a label for each cluster
#' @param spat_center_point_size size of the center point
#' @param spat_center_point_border_col border color of spatial center points
#' @param spat_center_point_border_stroke border strike size of spatial center
#' points
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param show_spatial_grid show spatial grid
#' @param spat_grid_name name of spatial grid to use
#' @param spat_grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param dim_other_point_size size of not selected dim cells
#' @param spat_other_point_size size of not selected spat cells
#' @param spat_other_cells_alpha alpha of not selected spat cells
#' @param dim_show_legend show legend of dimension reduction plot
#' @param spat_show_legend show legend of spatial plot
#' @param dim_background_color background color of points in dim. reduction
#' space
#' @param spat_background_color background color of spatial points
#' @param vor_border_color border color for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparency of voronoi 'cells'
#' @details Description of parameters.
#' @family spatial and dimension reduction visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatDimPlot2D(g)
#'
#' @export
#' @seealso \code{\link{spatDimPlot3D}}
spatDimPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        spat_loc_name = NULL,
        plot_alignment = c("vertical", "horizontal"),
        dim_reduction_to_use = "umap",
        dim_reduction_name = NULL,
        dim1_to_use = 1,
        dim2_to_use = 2,
        sdimx = "sdimx",
        sdimy = "sdimy",
        spat_enr_names = NULL,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        dim_point_shape = c("border", "no_border"),
        dim_point_size = 1,
        dim_point_alpha = 1,
        dim_point_border_col = "black",
        dim_point_border_stroke = 0.1,
        spat_point_shape = c("border", "no_border", "voronoi"),
        spat_point_size = 1,
        spat_point_alpha = 1,
        spat_point_border_col = "black",
        spat_point_border_stroke = 0.1,
        dim_show_cluster_center = FALSE,
        dim_show_center_label = TRUE,
        dim_center_point_size = 4,
        dim_center_point_border_col = "black",
        dim_center_point_border_stroke = 0.1,
        dim_label_size = 4,
        dim_label_fontface = "bold",
        spat_show_cluster_center = FALSE,
        spat_show_center_label = FALSE,
        spat_center_point_size = 4,
        spat_center_point_border_col = "blue",
        spat_center_point_border_stroke = 0.1,
        spat_label_size = 4,
        spat_label_fontface = "bold",
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        nn_network_alpha = 0.05,
        show_spatial_network = FALSE,
        spat_network_name = "Delaunay_network",
        spat_network_color = "blue",
        spat_network_alpha = 0.5,
        show_spatial_grid = FALSE,
        spat_grid_name = "spatial_grid",
        spat_grid_color = "blue",
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        dim_other_point_size = 1,
        spat_other_point_size = 1,
        spat_other_cells_alpha = 0.5,
        dim_show_legend = FALSE,
        spat_show_legend = FALSE,
        legend_text = 10,
        legend_symbol_size = 2,
        dim_background_color = "white",
        spat_background_color = "white",
        vor_border_color = "white",
        vor_max_radius = 200,
        vor_alpha = 1,
        axis_text = 8,
        axis_title = 8,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatDimPlot2D",
        view = NULL,
        space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "dimension_reduction",
            "images"))
    # deprecation message
    if (!is.null(largeImage_name)) {
        deprecate_warn(
            when = "0.2.0",
            what = "spatDimPlot2D(largeImage_name)",
            details = c(
                "Use `image_name` argument instead for all images to plot."
            )
        )
        image_name <- c(image_name, largeImage_name)
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

    plot_alignment <- match.arg(plot_alignment,
        choices = c("vertical", "horizontal")
    )


    # create matching cell_color_code
    if (is.null(cell_color_code)) {
        if (is.character(cell_color)) {
            cell_metadata <- pDataDT(gobject,
                spat_unit = spat_unit,
                feat_type = feat_type
            )
            if (cell_color %in% colnames(cell_metadata)) {
                if (color_as_factor == TRUE) {
                    number_colors <- length(
                        unique(cell_metadata[[cell_color]])
                    )
                    cell_color_code <- set_default_color_discrete_cell(
                        instrs = instructions(gobject)
                    )(n = number_colors)
                    names(cell_color_code) <- unique(
                        cell_metadata[[cell_color]]
                    )
                    cell_color_code <- cell_color_code
                }
            }
        }
    }

    # dimension reduction plot
    dmpl <- dimPlot2D(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        group_by = NULL,
        group_by_subset = NULL,
        dim_reduction_to_use = dim_reduction_to_use,
        dim_reduction_name = dim_reduction_name,
        dim1_to_use = dim1_to_use,
        dim2_to_use = dim2_to_use,
        spat_enr_names = spat_enr_names,
        cell_color = cell_color,
        color_as_factor = color_as_factor,
        cell_color_code = cell_color_code,
        cell_color_gradient = cell_color_gradient,
        gradient_midpoint = gradient_midpoint,
        gradient_style = gradient_style,
        gradient_limits = gradient_limits,
        select_cell_groups = select_cell_groups,
        select_cells = select_cells,
        point_shape = dim_point_shape,
        point_size = dim_point_size,
        point_alpha = dim_point_alpha,
        point_border_col = dim_point_border_col,
        point_border_stroke = dim_point_border_stroke,
        show_cluster_center = dim_show_cluster_center,
        show_center_label = dim_show_center_label,
        center_point_size = dim_center_point_size,
        center_point_border_col = dim_center_point_border_col,
        center_point_border_stroke = dim_center_point_border_stroke,
        label_size = dim_label_size,
        label_fontface = dim_label_fontface,
        show_NN_network = show_NN_network,
        nn_network_to_use = nn_network_to_use,
        network_name = network_name,
        edge_alpha = nn_network_alpha,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        other_point_size = dim_other_point_size,
        show_legend = dim_show_legend,
        legend_text = legend_text,
        legend_symbol_size = legend_symbol_size,
        background_color = dim_background_color,
        axis_text = axis_text,
        axis_title = axis_title,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE
    )

    # spatial plot
    spl <- spatPlot2D(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        show_image = show_image,
        gimage = gimage,
        image_name = image_name,
        spat_loc_name = spat_loc_name,
        group_by = NULL,
        group_by_subset = NULL,
        sdimx = sdimx,
        sdimy = sdimy,
        spat_enr_names = spat_enr_names,
        cell_color = cell_color,
        cell_color_code = cell_color_code,
        color_as_factor = color_as_factor,
        cell_color_gradient = cell_color_gradient,
        gradient_midpoint = gradient_midpoint,
        gradient_style = gradient_style,
        gradient_limits = gradient_limits,
        select_cell_groups = select_cell_groups,
        select_cells = select_cells,
        point_shape = spat_point_shape,
        point_size = spat_point_size,
        point_alpha = spat_point_alpha,
        point_border_col = spat_point_border_col,
        point_border_stroke = spat_point_border_stroke,
        show_cluster_center = spat_show_cluster_center,
        show_center_label = spat_show_center_label,
        center_point_size = spat_center_point_size,
        center_point_border_col = spat_center_point_border_col,
        center_point_border_stroke = spat_center_point_border_stroke,
        label_size = spat_label_size,
        label_fontface = spat_label_fontface,
        show_network = show_spatial_network,
        spatial_network_name = spat_network_name,
        network_color = spat_network_color,
        network_alpha = spat_network_alpha,
        show_grid = show_spatial_grid,
        spatial_grid_name = spat_grid_name,
        grid_color = spat_grid_color,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        other_point_size = spat_other_point_size,
        other_cells_alpha = spat_other_cells_alpha,
        coord_fix_ratio = 1,
        title = "",
        show_legend = spat_show_legend,
        legend_text = legend_text,
        legend_symbol_size = legend_symbol_size,
        background_color = spat_background_color,
        vor_border_color = vor_border_color,
        vor_max_radius = vor_max_radius,
        vor_alpha = vor_alpha,
        axis_text = axis_text,
        axis_title = axis_title,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE
    )


    if (plot_alignment == "vertical") {
        ncol <- 1
        nrow <- 2
        combo_plot <- cowplot::plot_grid(dmpl, spl,
            ncol = ncol,
            nrow = nrow, rel_heights = c(1),
            rel_widths = c(1), align = "v"
        )
    } else {
        ncol <- 2
        nrow <- 1
        combo_plot <- cowplot::plot_grid(dmpl, spl,
            ncol = ncol,
            nrow = nrow, rel_heights = c(1),
            rel_widths = c(1), align = "h"
        )
    }

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




#' @rdname spatDimPlot
#' @param \dots spatDimPlot(...) passes to spatDimPlot2D()
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatDimPlot(g)
#'
#' @export
spatDimPlot <- function(gobject, ...) {
    spatDimPlot2D(gobject, ...)
}



## ** spatial feature plotting ####

#' @title spatFeatPlot2D_single
#' @name spatFeatPlot2D_single
#' @description Visualize cells and feature expression according to
#' spatial coordinates
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_image_params
#' @inheritParams plot_params
#' @param plot_method character. How to plot the points. Either "ggplot" for
#' the default or "scattermore" for a faster rasterized option
#' @param largeImage_name deprecated
#' @param spat_loc_name name of spatial locations
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param expression_values gene expression values to use
#' @param feats features to show
#' @param order order points according to feature expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param edge_alpha alpha of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with
#' ggplot alpha parameter
#' @param point_shape shape of points (border, no_border or voronoi)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of points
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param coord_fix_ratio fix ratio between x and y-axis (default = 1)
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparency of voronoi 'cells'
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @details Description of parameters.
#' @family spatial feature expression visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatFeatPlot2D_single(g, feats = c("Gna12", "Ccnd2", "Btbd17"))
#'
#' @export
#' @seealso \code{\link{spatFeatPlot3D}}
spatFeatPlot2D_single <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        plot_method = "ggplot",
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        spat_loc_name = "raw",
        sdimx = "sdimx",
        sdimy = "sdimy",
        spat_enr_names = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        feats,
        order = TRUE,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        show_network = FALSE,
        network_color = NULL,
        edge_alpha = 0.5,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        midpoint = 0,
        scale_alpha_with_expression = FALSE,
        point_shape = c("border", "no_border", "voronoi"),
        point_size = 1,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        coord_fix_ratio = 1,
        show_legend = TRUE,
        show_axes = NULL,
        legend_text = 8,
        background_color = "white",
        vor_border_color = "white",
        vor_alpha = 1,
        vor_max_radius = 200,
        axis_text = 8,
        axis_title = 8,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        theme_param = list(),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatFeatPlot2D_single",
        view = NULL,
        space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "images"))
    # data.table variables
    cell_ID <- NULL

    # deprecation message
    if (!is.null(largeImage_name)) {
        deprecate_warn(
            when = "0.2.0",
            what = "spatFeatPlot2D_single(largeImage_name)",
            details = c(
                "Use `image_name` argument instead for all images to plot."
            )
        )
        image_name <- c(image_name, largeImage_name)
    }

    show_axes <- show_axes %null% TRUE

    # print, return and save parameters
    show_plot <- ifelse(is.null(show_plot),
        instructions(gobject, param = "show_plot"),
        show_plot
    )
    save_plot <- ifelse(is.null(save_plot),
        instructions(gobject, param = "save_plot"),
        save_plot
    )
    return_plot <- ifelse(is.null(return_plot),
        instructions(gobject, param = "return_plot"),
        return_plot
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

    ## get giotto image(s) ##
    if (isTRUE(show_image) && is.null(gimage)) {
        gimage <- getGiottoImage(
            gobject = gobject,
            name = image_name
        )
    }

    # point shape
    point_shape <- match.arg(point_shape,
        choices = c("border", "no_border", "voronoi")
    )

    # expression values
    values <- match.arg(
        expression_values,
        unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
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
    selected_feats <- setdiff(colnames(t_sub_expr_data_DT), "cell_ID")


    ## extract cell locations
    if (is.null(spat_loc_name)) {
        if (!is.null(slot(gobject, "spatial_locs"))) {
            spat_loc_name <- list_spatial_locations_names(
                gobject,
                spat_unit = spat_unit
            )[[1]]
        } else {
            spat_loc_name <- NULL
            warning("No spatial locations have been found")
            return(NULL)
        }
    }

    cell_locations <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        name = spat_loc_name,
        output = "data.table",
        copy_obj = TRUE
    )
    if (is.null(cell_locations)) {
        return(NULL)
    }

    ## extract spatial network
    if (show_network) {
        spatial_network <- getSpatialNetwork(
            gobject = gobject,
            spat_unit = spat_unit,
            name = spatial_network_name,
            output = "networkDT"
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

    cell_locations_metadata_feats <- merge(
        cell_locations_metadata,
        t_sub_expr_data_DT,
        by = "cell_ID"
    )


    ## plotting ##
    savelist <- list()

    for (feat in selected_feats) {
        # order spatial units (e.g. cell IDs) based on expression of feature
        if (isTRUE(order)) {
            cell_locations_metadata_feats <- cell_locations_metadata_feats[
                order(get(feat))
            ]
        }


        pl <- ggplot2::ggplot()
        pl <- pl + ggplot2::theme_classic()


        ## plot image ## TODO
        ## plot image ##
        if (isTRUE(show_image) && !is.null(gimage)) {
            pl <- plot_spat_image_layer_ggplot(
                gg_obj = pl,
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                spat_loc_name = spat_loc_name,
                gimage = gimage
            )
        }

        ## plot network or grid first if point_shape is border or no_border
        ## point
        if (point_shape %in% c("border", "no_border")) {
            ## plot spatial network
            if (!is.null(spatial_network) && isTRUE(show_network)) {
                edge_alpha <- edge_alpha %null% 0.5
                network_color <- network_color %null% "red"
                xbegin <- paste0(sdimx, "_begin")
                ybegin <- paste0(sdimy, "_begin")
                xend <- paste0(sdimx, "_end")
                yend <- paste0(sdimy, "_end")
                pl <- pl + ggplot2::geom_segment(
                    data = spatial_network,
                    aes_string2(
                        x = xbegin,
                        y = ybegin,
                        xend = xend,
                        yend = yend
                    ),
                    color = network_color,
                    linewidth = 0.5,
                    alpha = edge_alpha
                )
            }

            ## plot spatial grid
            if (!is.null(spatial_grid) && isTRUE(show_grid)) {
                if (is.null(grid_color)) grid_color <- "black"

                xmin <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimx
                ), "_start")
                ymin <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimy
                ), "_start")
                xmax <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimx
                ), "_end")
                ymax <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimy
                ), "_end")

                pl <- pl + ggplot2::geom_rect(
                    data = spatial_grid,
                    aes_string2(
                        xmin = xmin,
                        xmax = xmax,
                        ymin = ymin,
                        ymax = ymax
                    ),
                    color = grid_color,
                    fill = NA
                )
            }
        }



        ### plot cells ###

        ## set gradient limits if needed ##
        if (!is.null(gradient_limits) &&
            is.vector(gradient_limits) &&
            length(gradient_limits) == 2) {
            cell_locations_metadata_feats[[feat]] <-
                scales::oob_squish(cell_locations_metadata_feats[[feat]],
                    range = gradient_limits
                )
        }

        if (is.null(gradient_midpoint)) {
            gradient_midpoint <- stats::median(
                cell_locations_metadata_feats[[feat]]
            )
        }


        if (point_shape %in% c("border", "no_border")) {
            # assemble points plotting params
            # * aes - dynamic values found in the `data`
            # * args - static values to set

            points_aes <- aes_string2(x = sdimx, y = sdimy)

            points_args <- list()
            # common args
            points_args$size <- point_size
            points_args$show.legend <- show_legend
            points_args$plot_method <- plot_method

            if (isTRUE(scale_alpha_with_expression)) {
                points_aes$alpha <- as.name(feat)
            } else {
                points_args$alpha <- point_alpha
            }

            switch(point_shape,
                "border" = {
                    points_aes$fill <- as.name(feat)

                    points_args$shape <- 21
                    points_args$colour <- point_border_col
                    points_args$stroke <- point_border_stroke
                    scale_type <- "fill"
                },
                "no_border" = {
                    points_aes$colour <- as.name(feat)

                    points_args$shape <- 19
                    points_args$ext <- ext(gobject, prefer = "spatlocs")
                    scale_type <- "color"
                }
            )

            # other data to add
            points_args$data <- cell_locations_metadata_feats
            points_args$mapping <- points_aes

            # add points
            pl <- pl + do.call(giotto_point, args = points_args)

            ## scale and labs ##
            pl <- pl + ggplot2::scale_alpha_continuous(guide = "none")
            pl <- pl + set_default_color_continuous_cell(
                colors = cell_color_gradient,
                instrs = instructions(gobject),
                midpoint = gradient_midpoint,
                style = gradient_style,
                guide = ggplot2::guide_colorbar(title = ""),
                type = scale_type
            )
            pl <- pl + ggplot2::labs(
                x = "coord x",
                y = "coord y",
                title = feat
            )
        }


        ## voronoi ##
        if (point_shape == "voronoi") {
            if (scale_alpha_with_expression == TRUE) {
                pl <- pl + ggforce::geom_voronoi_tile(
                    data = cell_locations_metadata_feats,
                    aes_string2(
                        x = sdimx,
                        y = sdimy,
                        group = -1L,
                        fill = feat,
                        alpha = feat
                    ),
                    colour = vor_border_color,
                    max.radius = vor_max_radius,
                    show.legend = show_legend
                )
            } else {
                pl <- pl + ggforce::geom_voronoi_tile(
                    data = cell_locations_metadata_feats,
                    aes_string2(
                        x = sdimx,
                        y = sdimy,
                        group = -1L,
                        fill = feat
                    ),
                    colour = vor_border_color,
                    max.radius = vor_max_radius,
                    show.legend = show_legend,
                    alpha = vor_alpha
                )
            }


            ## plot spatial network
            if (!is.null(spatial_network) && show_network == TRUE) {
                if (is.null(network_color)) {
                    network_color <- "red"
                }
                xbegin <- paste0(sdimx, "_begin")
                ybegin <- paste0(sdimy, "_begin")
                xend <- paste0(sdimx, "_end")
                yend <- paste0(sdimy, "_end")
                pl <- pl + ggplot2::geom_segment(
                    data = spatial_network,
                    aes_string2(
                        x = xbegin,
                        y = ybegin,
                        xend = xend,
                        yend = yend
                    ),
                    color = network_color,
                    linewidth = 0.5,
                    alpha = 0.5
                )
            }

            ## plot spatial grid
            if (!is.null(spatial_grid) & show_grid == TRUE) {
                if (is.null(grid_color)) grid_color <- "black"

                xmin <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimx
                ), "_start")
                ymin <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimy
                ), "_start")
                xmax <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimx
                ), "_end")
                ymax <- paste0(gsub(
                    pattern = "sdim",
                    replacement = "", x = sdimy
                ), "_end")

                pl <- pl + ggplot2::geom_rect(
                    data = spatial_grid,
                    aes_string2(
                        xmin = xmin, xmax = xmax,
                        ymin = ymin, ymax = ymax
                    ),
                    color = grid_color,
                    fill = NA
                )
            }


            ## scale and labs ##
            pl <- pl + ggplot2::scale_alpha_continuous(guide = "none")
            pl <- pl + set_default_color_continuous_cell(
                colors = cell_color_gradient,
                instrs = instructions(gobject),
                midpoint = gradient_midpoint,
                style = gradient_style,
                guide = ggplot2::guide_colorbar(title = ""),
                type = "fill"
            )
            pl <- pl + ggplot2::labs(x = "coord x", y = "coord y", title = feat)
        }

        ## adjust theme setting
        gg_theme_args <- c(
            theme_param,
            legend_text = legend_text,
            axis_title = axis_title,
            axis_text = axis_text,
            background_color = background_color
        )
        pl <- pl + do.call(.gg_theme, args = gg_theme_args)

        if (!show_axes) {
            pl <- .theme_remove_axes(pl)
        }

        pl <- .aspect_ratio(pl, coord_fix_ratio)

        savelist[[feat]] <- pl
    }

    # combine plots with cowplot
    combo_plot <- cowplot::plot_grid(
        plotlist = savelist,
        ncol = set_default_cow_n_col(
            cow_n_col = cow_n_col,
            nr_plots = length(savelist)
        ),
        rel_heights = cow_rel_h,
        rel_widths = cow_rel_w,
        align = cow_align
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


#' @title Plot data in physical space 2D
#' @name spatFeatPlot2D
#' @description Visualize cells and feature expression according to
#' spatial coordinates
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_image_params
#' @inheritParams plot_params
#' @inheritParams plot_spatnet_params
#' @param plot_method method to plot points. Either "ggplot" (default) or
#' "scattermore" (rasterized and faster for large datasets)
#' @param largeImage_name deprecated
#' @param spat_loc_name name of spatial locations
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param expression_values gene expression values to use
#' @param feats features to show
#' @param order order points according to feature expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param edge_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha
#' parameter
#' @param coord_fix_ratio fix ratio between x and y-axis (default = 1)
#' @param background_color color of plot background
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparency of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @details Description of parameters.
#' @family spatial feature expression visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatFeatPlot2D(g, feats = "Gna12")
#'
#' @export
#' @seealso \code{\link{spatFeatPlot3D}}
spatFeatPlot2D <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        plot_method = "ggplot",
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        spat_loc_name = NULL,
        group_by = NULL,
        group_by_subset = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        expression_values = c("normalized", "scaled", "custom"),
        feats,
        order = TRUE,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        show_network = FALSE,
        network_color = NULL,
        edge_alpha = NULL,
        spatial_network_name = "Delaunay_network",
        show_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        midpoint = 0,
        scale_alpha_with_expression = FALSE,
        point_shape = c("border", "no_border", "voronoi"),
        point_size = 1,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        coord_fix_ratio = 1,
        show_legend = TRUE,
        show_axes = NULL,
        legend_text = 8,
        background_color = "white",
        vor_border_color = "white",
        vor_alpha = 1,
        vor_max_radius = 200,
        axis_text = 8,
        axis_title = 8,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        theme_param = list(),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatFeatPlot2D",
        view = NULL,
        space = NULL) {
    if (inherits(gobject, "giottoMulti")) {
        return(.gg_multi_dispatch_spatial(
            plot_fn = spatFeatPlot2D,
            named = mget(names(formals())),
            gobject = gobject, view = view, space = space
        ))
    }
    .gg_assert_giotto_single(gobject)
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "images"))
    # deprecation message
    if (!is.null(largeImage_name)) {
        deprecate_warn(
            when = "0.2.0",
            what = "spatFeatPlot2D(largeImage_name)",
            details = c(
                "Use `image_name` argument instead for all images to plot."
            )
        )
        image_name <- c(image_name, largeImage_name)
    }

    plot_method <- match.arg(plot_method, c("ggplot", "scattermore"))

    point_shape <- match.arg(point_shape, c("border", "no_border", "voronoi"))

    if (identical(plot_method, "scattermore") && point_shape != "no_border") {
        warning("point_shape changed to \"no_border\" for scattermore")
        point_shape <- "no_border"
    }

    # create args list needed for each call to spatFeatPlot2D_single()
    # 1. - grab all params available
    # 2. - subset to those needed
    sfp_params <- get_args_list()
    sfp_params <- sfp_params[c(
        # [gobject params]
        "gobject", "feat_type", "spat_unit",
        # [image params]
        "show_image", "gimage", "image_name",
        # [spatlocs params]
        "spat_loc_name", "sdimx", "sdimy",
        # [expression params]
        "expression_values", "feats", "order",
        # [point aes]
        "cell_color_gradient", "gradient_midpoint", "gradient_style",
        "gradient_limits", "midpoint", "scale_alpha_with_expression",
        "point_shape", "plot_method",
        "point_size", "point_alpha", "point_border_col", "point_border_stroke",
        # [voronoi-point params]
        "vor_border_color", "vor_alpha", "vor_max_radius",
        # [network aes]
        "show_network", "network_color", "edge_alpha", "spatial_network_name",
        # [grid aes]
        "show_grid", "grid_color", "spatial_grid_name",
        # [figure params]
        "coord_fix_ratio", "show_legend", "legend_text", "background_color",
        "axis_text", "axis_title", "show_axes",
        "cow_n_col", "cow_rel_h", "cow_rel_w", "cow_align",
        # [return params]
        "show_plot", "return_plot", "save_plot", "save_param",
        "default_save_name",
        # [theme params]
        "theme_param"
    )]

    ## check group_by
    if (is.null(group_by)) { # ----------------------------------------------- #

        do.call(spatFeatPlot2D_single, args = sfp_params)
    } else { # -------------------------------------------------------------- #

        # Set feat_type and spat_unit
        spat_unit <- set_default_spat_unit(
            gobject = gobject,
            spat_unit = spat_unit
        )
        feat_type <- set_default_feat_type(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit
        )
        # ! update spat_unit & feat_type in static params ! #

        ## check metadata for valid group_by information
        comb_metadata <- combineMetadata(
            gobject = gobject,
            spat_unit = spat_unit,
            spat_loc_name = spat_loc_name,
            feat_type = feat_type
        )
        possible_meta_groups <- colnames(comb_metadata)

        ## error if group_by col is not found
        if (!group_by %in% possible_meta_groups) {
            stop("group_by ", group_by, " was not found in pDataDT()")
        }

        unique_groups <- unique(comb_metadata[[group_by]])

        # subset unique_groups
        # These unique_groups will be used to iterate through subsetting then
        # plotting the giotto object multiple times.
        if (!is.null(group_by_subset)) {
            not_found <- group_by_subset[!group_by_subset %in% unique_groups]
            if (length(not_found) > 0) {
                message("the following subset was not found: ", not_found)
            }
            unique_groups <- unique_groups[unique_groups %in% group_by_subset]
        }


        # group_by static settings #
        # update these params
        sfp_params$spat_unit <- spat_unit
        sfp_params$feat_type <- feat_type
        # apply group_by specific defaults
        sfp_params$cow_n_col <- 1
        sfp_params$show_plot <- FALSE
        sfp_params$return_plot <- TRUE
        sfp_params$save_plot <- FALSE
        sfp_params$default_save_name <- "spatFeatPlot2D"


        ## plotting ##
        savelist <- list()

        for (group_id in seq_along(unique_groups)) {
            group <- unique_groups[group_id]

            subset_cell_IDs <- comb_metadata[
                get(group_by) == group
            ][["cell_ID"]]
            sfp_params$gobject <- subsetGiotto(
                gobject = gobject,
                feat_type = feat_type,
                spat_unit = spat_unit,
                cell_ids = subset_cell_IDs,
                verbose = FALSE
            )

            # use a different image per group if there are the same number of
            # names provided as there are groups
            # Otherwise, use the same image (or NULL) for all groups (default)
            if (length(unique_groups) == length(image_name)) {
                sfp_params$image_name <- image_name[group_id]
            }


            pl <- do.call(spatFeatPlot2D_single, args = sfp_params)

            savelist[[group_id]] <- pl
        }

        # combine plots with cowplot
        combo_plot <- cowplot::plot_grid(
            plotlist = savelist,
            ncol = set_default_cow_n_col(
                cow_n_col = cow_n_col,
                nr_plots = length(savelist)
            ),
            rel_heights = cow_rel_h,
            rel_widths = cow_rel_w,
            align = cow_align
        )

        # output
        return(
            plot_output_handler(
                gobject = gobject,
                plot_object = combo_plot,
                save_plot = save_plot,
                return_plot = return_plot,
                show_plot = show_plot,
                default_save_name = default_save_name,
                save_param = save_param,
                else_return = NULL
            )
        )
    } # --------------------------------------------------------------------- #
}








## ** dim reduction feature plotting ####

.dimFeatPlot2D_single <- function(selected_feats, ...) {
    lapply(selected_feats, function(feat) {
        .dimFeatPlot_single_feat(feat, ...)
    })
}

.dimFeatPlot_single_feat <- function(
        data,
        feat,
        feats,
        dim_names,
        order = TRUE,
        group_id = NULL,
        show_NN_network = FALSE,
        network_color = NULL,
        from_dim_names = NULL,
        to_dim_names = NULL,
        annotated_network_DT = NULL,
        edge_alpha = NULL,
        scale_alpha_with_expression = FALSE,
        point_shape,
        point_size = 1,
        point_alpha = 1,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        point_border_col = "black",
        point_border_stroke = 0.1,
        show_legend = TRUE,
        legend_text = 10,
        background_color = "white",
        axis_text = 8,
        axis_title = 8,
        instrs) {

    to <- from <- NULL # NSE vars

    # order spatial units (e.g. cell IDs) based on expression of feature
    if (isTRUE(order)) {
        data <- data[order(get(feat))]
    }

    ## OLD need to be combined ##
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    # network layer
    if (show_NN_network == TRUE) {

        annotated_network_DT <- annotated_network_DT[
            to %in% data$cell_ID & from %in% data$cell_ID
        ]


        if (is.null(edge_alpha)) {
            edge_alpha <- 0.5
            pl <- pl + ggplot2::geom_segment(
                data = annotated_network_DT,
                aes_string2(
                    x = from_dim_names[1], y = from_dim_names[2],
                    xend = to_dim_names[1], yend = to_dim_names[2]
                ),
                alpha = edge_alpha,
                color = network_color,
                linewidth = 0.1,
                show.legend = FALSE
            )
        } else if (is.numeric(edge_alpha)) {
            pl <- pl + ggplot2::geom_segment(
                data = annotated_network_DT,
                aes_string2(
                    x = from_dim_names[1], y = from_dim_names[2],
                    xend = to_dim_names[1], yend = to_dim_names[2]
                ),
                alpha = edge_alpha,
                color = network_color,
                linewidth = 0.1,
                show.legend = FALSE
            )
        } else if (is.character(edge_alpha)) {
            if (edge_alpha %in% colnames(annotated_network_DT)) {
                pl <- pl + ggplot2::geom_segment(
                    data = annotated_network_DT,
                    aes_string2(
                        x = from_dim_names[1],
                        y = from_dim_names[2],
                        xend = to_dim_names[1],
                        yend = to_dim_names[2],
                        alpha = edge_alpha
                    ),
                    color = network_color,
                    show.legend = FALSE
                )
            }
        }
    }


    ## point layer ##
    if (is.null(feats)) {
        cell_color <- "lightblue"
        message("no feats selected")
        pl <- pl + ggplot2::geom_point(
            data = data,
            aes_string2(x = dim_names[1], dim_names[2]),
            fill = cell_color,
            show.legend = show_legend,
            size = point_size,
            alpha = point_alpha
        )
    } else {
        ## set gradient limits if needed ##
        if (!is.null(gradient_limits) & is.vector(gradient_limits) &
            length(gradient_limits) == 2) {
            lower_lim <- gradient_limits[[1]]
            upper_lim <- gradient_limits[[2]]
            numeric_data <- data[[feat]]
            limit_numeric_data <- ifelse(
                numeric_data > upper_lim,
                upper_lim,
                ifelse(numeric_data < lower_lim, lower_lim, numeric_data)
            )
            data[[feat]] <- limit_numeric_data
        }

        if (is.null(gradient_midpoint)) {
            gradient_midpoint <- stats::median(data[[feat]])
        }



        ## with border ##
        if (point_shape == "border") {
            if (scale_alpha_with_expression == TRUE) {
                pl <- pl + ggplot2::geom_point(
                    data = data,
                    aes_string2(
                        x = dim_names[1],
                        y = dim_names[2],
                        fill = feat,
                        alpha = feat
                    ),
                    show.legend = show_legend,
                    shape = 21,
                    size = point_size,
                    color = point_border_col,
                    stroke = point_border_stroke
                )
            } else {
                pl <- pl + ggplot2::geom_point(
                    data = data,
                    aes_string2(
                        x = dim_names[1],
                        y = dim_names[2],
                        fill = feat
                    ),
                    show.legend = show_legend,
                    shape = 21,
                    size = point_size,
                    color = point_border_col,
                    stroke = point_border_stroke,
                    alpha = point_alpha
                )
            }

            ## scale and labs ##
            pl <- pl + ggplot2::scale_alpha_continuous(guide = "none")
            pl <- pl + set_default_color_continuous_cell(
                colors = cell_color_gradient,
                instrs = instrs,
                midpoint = gradient_midpoint,
                style = gradient_style,
                guide = ggplot2::guide_colorbar(title = ""),
                type = "fill"
            )
        }


        ## without border ##
        if (point_shape == "no_border") {
            if (scale_alpha_with_expression == TRUE) {
                pl <- pl + ggplot2::geom_point(
                    data = data,
                    aes_string2(
                        x = dim_names[1],
                        y = dim_names[2],
                        color = feat,
                        alpha = feat
                    ),
                    show.legend = show_legend,
                    shape = 19,
                    size = point_size
                )
            } else {
                pl <- pl + ggplot2::geom_point(
                    data = data,
                    aes_string2(
                        x = dim_names[1],
                        y = dim_names[2],
                        color = feat
                    ),
                    show.legend = show_legend,
                    shape = 19,
                    size = point_size,
                    alpha = point_alpha
                )
            }

            ## scale and labs ##
            pl <- pl + ggplot2::scale_alpha_continuous(guide = "none")
            pl <- pl + set_default_color_continuous_cell(
                colors = cell_color_gradient,
                instrs = instrs,
                midpoint = gradient_midpoint,
                style = gradient_style,
                guide = ggplot2::guide_colorbar(title = ""),
                type = "color"
            )
        }
    }

    ## add title
    pl <- pl + ggplot2::labs(
        x = "coord x",
        y = "coord y",
        title = paste(c(group_id, feat), collapse = "-")
    )

    ## aesthetics
    pl <- pl + ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_text),
        axis.title = element_text(size = axis_title),
        axis.text = element_text(size = axis_text),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = background_color)
    )

    return(pl)
}



#' @title dimFeatPlot2D
#' @name dimFeatPlot2D
#' @description Visualize gene expression according to dimension reduction
#' coordinates
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_nn_net_params
#' @inheritParams plot_dimred_params
#' @inheritParams plot_params
#' @param expression_values gene expression values to use
#' @param feats features to show
#' @param order order points according to feature expression
#' @param scale_alpha_with_expression scale expression with ggplot alpha
#' parameter
#' @details Description of parameters.
#' @family dimension reduction feature expression visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' dimFeatPlot2D(g, feats = c("Gna12", "Ccnd2", "Btbd17"))
#'
#' # with group_by
#' dimFeatPlot2D(g,
#'     feats = c("Gna12"),
#'     group_by = "leiden_clus",
#'     gradient_midpoint = 3 # setting a specific midpoint can be helpful
#' )
#' # with group_by and group_by_subset
#' dimFeatPlot2D(g,
#'     feats = c("Gna12", "Ccnd2", "Btbd17"),
#'     group_by = "leiden_clus",
#'     group_by_subset = c(2, 5)
#' )
#' @export
dimFeatPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        feats = NULL,
        order = TRUE,
        group_by = NULL,
        group_by_subset = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = NULL,
        dim1_to_use = 1,
        dim2_to_use = 2,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        network_color = "lightgray",
        edge_alpha = NULL,
        scale_alpha_with_expression = FALSE,
        point_shape = c("border", "no_border"),
        point_size = 1,
        point_alpha = 1,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        point_border_col = "black",
        point_border_stroke = 0.1,
        show_legend = TRUE,
        legend_text = 10,
        background_color = "white",
        axis_text = 8,
        axis_title = 8,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "dimFeatPlot2D",
        view = NULL,
        space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "dimension_reduction",
            "spatial_enrichment", "expression"))

    handle_errors({
        # print, return and save parameters
        show_plot <- ifelse(is.null(show_plot),
            instructions(gobject, param = "show_plot"),
            show_plot
        )
        save_plot <- ifelse(is.null(save_plot),
            instructions(gobject, param = "save_plot"),
            save_plot
        )
        return_plot <- ifelse(is.null(return_plot),
            instructions(gobject, param = "return_plot"),
            return_plot
        )

        # point shape
        point_shape <- match.arg(
            point_shape, choices = c("border", "no_border"))

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

        # specify dim_reduction_name according to provided feat_type
        if (!is.null(dim_reduction_to_use)) {
            if (is.null(dim_reduction_name)) {
                if (feat_type == "rna") {
                    dim_reduction_name <- dim_reduction_to_use
                } else {
                    dim_reduction_name <- paste0(
                        feat_type, ".",
                        dim_reduction_to_use
                    )
                }
            }
        }


        # expression values
        values <- match.arg(
            expression_values,
            unique(c(
                "normalized", "scaled", "custom",
                expression_values
            ))
        )

        # data.table variables
        cell_ID <- NULL

        if (length(feats) == 0) {
            stop("No `feats` selected to plot.", call. = FALSE)
        }

        # collect feature values using spatValues
        t_sub_expr_data_DT <- spatValues(
            gobject,
            feats = feats,
            spat_unit = spat_unit,
            feat_type = feat_type,
            expression_values = values,
            verbose = FALSE
        )
        selected_feats <- setdiff(colnames(t_sub_expr_data_DT), "cell_ID")
        if (length(selected_feats) == 0) {
            stop("Selected `feats` not found in expression information",
                call. = FALSE
            )
        }


        ## dimension reduction ##
        dim_dfr <- getDimReduction(
            gobject = gobject,
            feat_type = feat_type,
            spat_unit = spat_unit,
            reduction = "cells",
            reduction_method = dim_reduction_to_use,
            name = dim_reduction_name,
            output = "data.table"
        )

        dim_names <- colnames(dim_dfr)
        dim_DT <- data.table::as.data.table(dim_dfr)
        dim_DT[, cell_ID := rownames(dim_dfr)]

        ## annotated cell metadata
        cell_metadata <- pDataDT(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type
        )

        annotated_DT <- data.table::merge.data.table(
            cell_metadata,
            dim_DT,
            by = "cell_ID"
        )

        ## merge feat info
        annotated_feat_DT <- data.table::merge.data.table(
            annotated_DT,
            t_sub_expr_data_DT,
            by = "cell_ID"
        )

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

            network_DT <- data.table::as.data.table(
                igraph::as_data_frame(selected_nn_network, what = "edges")
            )

            # annotated network
            old_dim_names <- dim_names

            annotated_network_DT <- data.table::merge.data.table(
                network_DT, dim_DT,
                by.x = "from", by.y = "cell_ID"
            )
            from_dim_names <- paste0("from_", old_dim_names)
            data.table::setnames(
                annotated_network_DT,
                old = old_dim_names,
                new = from_dim_names
            )

            annotated_network_DT <- data.table::merge.data.table(
                annotated_network_DT, dim_DT,
                by.x = "to", by.y = "cell_ID"
            )
            to_dim_names <- paste0("to_", old_dim_names)
            data.table::setnames(
                annotated_network_DT,
                old = old_dim_names,
                new = to_dim_names
            )
        }

        # params list for plotting
        a <- list(
            feats = feats,
            selected_feats = selected_feats,
            order = order,
            dim_names = dim_names,
            show_NN_network = show_NN_network,
            edge_alpha = edge_alpha,
            scale_alpha_with_expression = scale_alpha_with_expression,
            point_shape = point_shape,
            point_size = point_size,
            point_alpha = point_alpha,
            cell_color_gradient = cell_color_gradient,
            gradient_midpoint = gradient_midpoint,
            gradient_style = gradient_style,
            gradient_limits = gradient_limits,
            point_border_col = point_border_col,
            point_border_stroke = point_border_stroke,
            show_legend = show_legend,
            legend_text = legend_text,
            background_color = background_color,
            axis_text = axis_text,
            axis_title = axis_title,
            instrs = instructions(gobject)
        )

        if (isTRUE(show_NN_network)) {
            a$network_color <- network_color
            a$from_dim_names <- from_dim_names
            a$to_dim_names <- to_dim_names
            a$annotated_network_DT <- annotated_network_DT
        }

        ## generate plot(s) ##
        if (is.null(group_by)) {
            a$data <- annotated_feat_DT
            savelist <- do.call(.dimFeatPlot2D_single, a)
        } else {
            datalist <- .groupby(annotated_feat_DT,
                group_by = group_by,
                group_by_subset = group_by_subset
            )
            gp_savelist <- lapply(datalist, function(data) {
                a$data <- data
                group_id <- data[[group_by]][[1]]
                a$group_id <- group_id
                do.call(.dimFeatPlot2D_single, a)
            })
            savelist <- do.call(c, gp_savelist)
        }


        # combine plots with cowplot
        combo_plot <- cowplot::plot_grid(
            plotlist = savelist,
            ncol = set_default_cow_n_col(
                cow_n_col = cow_n_col,
                nr_plots = length(savelist)
            ),
            rel_heights = cow_rel_h,
            rel_widths = cow_rel_w,
            align = cow_align
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
                    gobject = gobject,
                    plot_object = combo_plot,
                    default_save_name = default_save_name
                ), save_param)
            )
        }

        ## return plot
        if (return_plot == TRUE) {
            return(combo_plot)
        }
    }, prefix = "dimFeatPlot2D()")
}






## ** spatial and dim reduction feature plotting ####


#' @title spatDimFeatPlot2D
#' @name spatDimFeatPlot2D
#' @description Visualize cells according to spatial AND dimension reduction
#' coordinates in ggplot mode
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_image_params
#' @inheritParams plot_params
#' @inheritParams plot_dimred_params
#' @inheritParams plot_spatnet_params
#' @inheritParams plot_nn_net_params
#' @param expression_values feat expression values to use
#' @param plot_alignment direction to align plot
#' @param feats features to show
#' @param order order points according to feature expression
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param dim_network_color color of NN network
#' @param dim_edge_alpha dim reduction plot: column to use for alpha of the
#' edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha
#' parameter
#' @param sdimx spatial x-axis dimension name (default = 'sdimx')
#' @param sdimy spatial y-axis dimension name (default = 'sdimy')
#' @param show_spatial_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spat_point_shape spatial points with border or
#' not (border or no_border)
#' @param spat_point_size spatial plot: point size
#' @param spat_point_alpha transparency of spatial points
#' @param spat_point_border_col color of border around points
#' @param spat_point_border_stroke stroke size of border around points
#' @param spat_edge_alpha edge alpha
#' @param dim_background_color color of plot background for dimension plot
#' @param spat_background_color color of plot background for spatial plot
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @details Description of parameters.
#' @family spatial and dimension reduction feature expression visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatDimFeatPlot2D(g, feats = c("Gna12", "Ccnd2", "Btbd17"))
#'
#' @export
spatDimFeatPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        expression_values = c("normalized", "scaled", "custom"),
        plot_alignment = c("vertical", "horizontal"),
        feats,
        order = TRUE,
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        dim_point_shape = c("border", "no_border"),
        dim_point_size = 1,
        dim_point_alpha = 1,
        dim_point_border_col = "black",
        dim_point_border_stroke = 0.1,
        show_NN_network = FALSE,
        show_spatial_network = FALSE,
        dim_network_color = "gray",
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        dim_edge_alpha = NULL,
        scale_alpha_with_expression = FALSE,
        sdimx = "sdimx",
        sdimy = "sdimy",
        spatial_network_name = "Delaunay_network",
        spatial_network_color = NULL,
        show_spatial_grid = FALSE,
        grid_color = NULL,
        spatial_grid_name = "spatial_grid",
        spat_point_shape = c("border", "no_border", "voronoi"),
        spat_point_size = 1,
        spat_point_alpha = 1,
        spat_point_border_col = "black",
        spat_point_border_stroke = 0.1,
        spat_edge_alpha = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        show_legend = TRUE,
        legend_text = 10,
        dim_background_color = "white",
        spat_background_color = "white",
        vor_border_color = "white",
        vor_max_radius = 200,
        vor_alpha = 1,
        axis_text = 8,
        axis_title = 8,
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatDimFeatPlot2D",
        view = NULL,
        space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "dimension_reduction",
            "images"))
    plot_alignment <- match.arg(plot_alignment,
        choices = c("vertical", "horizontal")
    )

    # dimension reduction plot
    dmpl <- dimFeatPlot2D(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        expression_values = expression_values,
        feats = feats,
        order = order,
        dim_reduction_to_use = dim_reduction_to_use,
        dim_reduction_name = dim_reduction_name,
        dim1_to_use = dim1_to_use,
        dim2_to_use = dim2_to_use,
        show_NN_network = show_NN_network,
        nn_network_to_use = nn_network_to_use,
        network_name = network_name,
        network_color = dim_network_color,
        edge_alpha = dim_edge_alpha,
        scale_alpha_with_expression = scale_alpha_with_expression,
        point_shape = dim_point_shape,
        point_size = dim_point_size,
        point_alpha = dim_point_alpha,
        cell_color_gradient = cell_color_gradient,
        gradient_midpoint = gradient_midpoint,
        gradient_style = gradient_style,
        gradient_limits = gradient_limits,
        point_border_col = dim_point_border_col,
        point_border_stroke = dim_point_border_stroke,
        show_legend = show_legend,
        legend_text = legend_text,
        background_color = dim_background_color,
        axis_text = axis_text,
        axis_title = axis_title,
        cow_n_col = cow_n_col,
        cow_rel_h = cow_rel_h,
        cow_rel_w = cow_rel_w,
        cow_align = cow_align,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE
    )

    # spatial plot
    spl <- spatFeatPlot2D(
        gobject = gobject,
        feat_type = feat_type,
        spat_unit = spat_unit,
        show_image = show_image,
        gimage = gimage,
        image_name = image_name,
        largeImage_name = largeImage_name,
        sdimx = sdimx,
        sdimy = sdimy,
        expression_values = expression_values,
        feats = feats,
        order = order,
        cell_color_gradient = cell_color_gradient,
        gradient_midpoint = gradient_midpoint,
        gradient_style = gradient_style,
        gradient_limits = gradient_limits,
        show_network = show_spatial_network,
        network_color = spatial_network_color,
        spatial_network_name = spatial_network_name,
        edge_alpha = spat_edge_alpha,
        show_grid = show_spatial_grid,
        grid_color = grid_color,
        spatial_grid_name = spatial_grid_name,
        scale_alpha_with_expression = scale_alpha_with_expression,
        point_shape = spat_point_shape,
        point_size = spat_point_size,
        point_alpha = spat_point_alpha,
        point_border_col = spat_point_border_col,
        point_border_stroke = spat_point_border_stroke,
        show_legend = show_legend,
        legend_text = legend_text,
        background_color = spat_background_color,
        vor_border_color = vor_border_color,
        vor_max_radius = vor_max_radius,
        vor_alpha = vor_alpha,
        axis_text = axis_text,
        axis_title = axis_title,
        cow_n_col = cow_n_col,
        cow_rel_h = cow_rel_h,
        cow_rel_w = cow_rel_w,
        cow_align = cow_align,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE
    )


    if (plot_alignment == "vertical") {
        ncol <- 1
        nrow <- 2
        combo_plot <- cowplot::plot_grid(
            dmpl,
            spl,
            ncol = ncol,
            nrow = nrow,
            rel_heights = c(1),
            rel_widths = c(1),
            align = "v"
        )
    } else {
        ncol <- 2
        nrow <- 1
        combo_plot <- cowplot::plot_grid(
            dmpl, spl,
            ncol = ncol,
            nrow = nrow,
            rel_heights = c(1),
            rel_widths = c(1),
            align = "h"
        )
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = combo_plot,
        show_plot = show_plot,
        save_plot = save_plot,
        return_plot = return_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}














#' @title spatCellPlot
#' @name spatCellPlot
#' @description Visualize cells according to spatial coordinates
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_params
#' @inheritParams plot_spatenr_params
#' @inheritParams plot_image_params
#' @inheritParams plot_spatnet_params
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_annotation_values numeric cell annotation columns
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparency of voronoi 'cells'
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @details Description of parameters.
#' @family spatial cell annotation visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatCellPlot2D(g, cell_annotation_values = "leiden_clus")
#'
#' @export
spatCellPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        sdimx = "sdimx",
        sdimy = "sdimy",
        spat_enr_names = NULL,
        cell_annotation_values = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        point_shape = c("border", "no_border", "voronoi"),
        point_size = 3,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        show_cluster_center = FALSE,
        show_center_label = FALSE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        show_network = FALSE,
        spatial_network_name = "Delaunay_network",
        network_color = NULL,
        network_alpha = 1,
        show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        grid_color = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        other_cells_alpha = 0.1,
        coord_fix_ratio = 1,
        show_legend = TRUE,
        legend_text = 8,
        legend_symbol_size = 1,
        background_color = "white",
        vor_border_color = "white",
        vor_max_radius = 200,
        vor_alpha = 1,
        axis_text = 8,
        axis_title = 8,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        theme_param = list(),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatCellPlot2D",
        view = NULL,
        space = NULL) {
    if (inherits(gobject, "giottoMulti")) {
        return(.gg_multi_dispatch_spatial(
            plot_fn = spatCellPlot2D,
            named = mget(names(formals())),
            gobject = gobject, view = view, space = space
        ))
    }
    .gg_assert_giotto_single(gobject)
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "images"))
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

    comb_metadata <- combineMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spat_enr_names = spat_enr_names
    )

    # keep only available columns
    possible_value_cols <- colnames(comb_metadata)
    if (is.null(cell_annotation_values)) {
        stop("you need to choose which continuous/numerical cell
            annotations or enrichments you want to visualize")
    }
    cell_annotation_values <- as.character(cell_annotation_values)
    cell_annotation_values <- cell_annotation_values[
        cell_annotation_values %in% possible_value_cols
    ]

    ## plotting ##
    savelist <- list()

    for (annot in cell_annotation_values) {
        pl <- spatPlot2D(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            show_image = show_image,
            gimage = gimage,
            image_name = image_name,
            largeImage_name = largeImage_name,
            group_by = NULL,
            group_by_subset = NULL,
            sdimx = sdimx,
            sdimy = sdimy,
            spat_enr_names = spat_enr_names,
            cell_color = annot,
            color_as_factor = FALSE,
            cell_color_gradient = cell_color_gradient,
            gradient_midpoint = gradient_midpoint,
            gradient_style = gradient_style,
            gradient_limits = gradient_limits,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            point_shape = point_shape,
            point_size = point_size,
            point_alpha = point_alpha,
            point_border_col = point_border_col,
            point_border_stroke = point_border_stroke,
            show_cluster_center = show_cluster_center,
            show_center_label = show_center_label,
            center_point_size = center_point_size,
            center_point_border_col = center_point_border_col,
            center_point_border_stroke = center_point_border_stroke,
            label_size = label_size,
            label_fontface = label_fontface,
            show_network = show_network,
            spatial_network_name = spatial_network_name,
            network_color = network_color,
            network_alpha = network_alpha,
            show_grid = show_grid,
            spatial_grid_name = spatial_grid_name,
            grid_color = grid_color,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            other_cells_alpha = other_cells_alpha,
            coord_fix_ratio = coord_fix_ratio,
            title = annot,
            show_legend = show_legend,
            legend_text = legend_text,
            legend_symbol_size = legend_symbol_size,
            background_color = background_color,
            vor_border_color = vor_border_color,
            vor_max_radius = vor_max_radius,
            vor_alpha = vor_alpha,
            axis_text = axis_text,
            axis_title = axis_title,
            theme_param = theme_param,
            # hardcoded on purpose below
            show_plot = FALSE,
            return_plot = TRUE,
            save_plot = FALSE,
            save_param = list(),
            default_save_name = "spatPlot2D"
        )


        savelist[[annot]] <- pl
    }


    # combine plots with cowplot
    combo_plot <- cowplot::plot_grid(
        plotlist = savelist,
        ncol = set_default_cow_n_col(
            cow_n_col = cow_n_col,
            nr_plots = length(savelist)
        ),
        rel_heights = cow_rel_h,
        rel_widths = cow_rel_w,
        align = cow_align
    )

    return(plot_output_handler(
        gobject = gobject,
        plot_object = combo_plot,
        save_plot = save_plot,
        show_plot = show_plot,
        return_plot = return_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}


#' @rdname spatCellPlot
#' @param \dots spatCellPlot(...) passes to spatCellPlot2D()
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatCellPlot(g, cell_annotation_values = "leiden_clus")
#'
#' @export
spatCellPlot <- function(...) {
    spatCellPlot2D(...)
}





#' @title dimCellPlot
#' @name dimCellPlot
#' @description Visualize cells according to dimension reduction coordinates.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_dimred_params
#' @inheritParams plot_nn_net_params
#' @inheritParams plot_spatenr_params
#' @inheritParams plot_params
#' @param cell_annotation_values numeric cell annotation columns
#' @details Description of parameters. For 3D plots see \code{\link{dimPlot3D}}
#' @family dimension reduction cell annotation visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' dimCellPlot2D(
#'     g,
#'     spat_enr_names = "cluster_metagene",
#'     cell_annotation_values = as.character(seq(4))
#' )
#'
#' @export
dimCellPlot2D <- function(
        gobject,
        spat_unit = NULL,
        feat_type = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        spat_enr_names = NULL,
        cell_annotation_values = NULL,
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        network_name = "sNN.pca",
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 0.5,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        edge_alpha = NULL,
        point_shape = c("border", "no_border"),
        point_size = 1,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        show_legend = TRUE,
        legend_text = 8,
        legend_symbol_size = 1,
        background_color = "white",
        axis_text = 8,
        axis_title = 8,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "dimCellPlot2D",
        view = NULL,
        space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "dimension_reduction",
            "spatial_enrichment", "expression"))
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

    comb_metadata <- combineMetadata(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spat_enr_names = spat_enr_names
    )

    # keep only available columns
    possible_value_cols <- colnames(comb_metadata)
    if (is.null(cell_annotation_values)) {
        stop("you need to choose which continuous/numerical cell annotations
            or enrichments you want to visualize")
    }
    cell_annotation_values <- cell_annotation_values[
        cell_annotation_values %in% possible_value_cols
    ]

    ## plotting ##
    savelist <- list()

    for (annot in cell_annotation_values) {
        pl <- dimPlot2D(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            group_by = NULL,
            group_by_subset = NULL,
            dim_reduction_to_use = dim_reduction_to_use,
            dim_reduction_name = dim_reduction_name,
            dim1_to_use = dim1_to_use,
            dim2_to_use = dim2_to_use,
            spat_enr_names = spat_enr_names,
            show_NN_network = show_NN_network,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            cell_color = annot,
            color_as_factor = FALSE,
            cell_color_code = cell_color_code,
            cell_color_gradient = cell_color_gradient,
            gradient_midpoint = gradient_midpoint,
            gradient_style = gradient_style,
            gradient_limits = gradient_limits,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_cluster_center = show_cluster_center,
            show_center_label = show_center_label,
            center_point_size = center_point_size,
            center_point_border_col = center_point_border_col,
            center_point_border_stroke = center_point_border_stroke,
            label_size = label_size,
            label_fontface = label_fontface,
            edge_alpha = edge_alpha,
            point_shape = point_shape,
            point_size = point_size,
            point_alpha = point_alpha,
            point_border_col = point_border_col,
            point_border_stroke = point_border_stroke,
            title = annot,
            show_legend = show_legend,
            legend_text = legend_text,
            legend_symbol_size = legend_symbol_size,
            background_color = background_color,
            axis_text = axis_text,
            axis_title = axis_title,
            show_plot = FALSE,
            return_plot = TRUE,
            save_plot = FALSE,
            save_param = list(),
            default_save_name = "dimPlot2D"
        )


        savelist[[annot]] <- pl
    }


    # combine plots with cowplot
    combo_plot <- cowplot::plot_grid(
        plotlist = savelist,
        ncol = set_default_cow_n_col(
            cow_n_col = cow_n_col,
            nr_plots = length(savelist)
        ),
        rel_heights = cow_rel_h,
        rel_widths = cow_rel_w,
        align = cow_align
    )

    return(plot_output_handler(
        gobject = gobject,
        plot_object = combo_plot,
        save_plot = save_plot,
        show_plot = show_plot,
        return_plot = return_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}




#' @rdname dimCellPlot
#' @param ... dimCellPlot(...) passes to dimCellPlot2D()
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' dimCellPlot(g, cell_annotation_values = "leiden_clus")
#'
#' @export
dimCellPlot <- function(gobject, ...) {
    dimCellPlot2D(gobject = gobject, ...)
}




#' @title spatDimCellPlot2D
#' @name spatDimCellPlot2D
#' @description Visualize numerical features of cells according to spatial
#' AND dimension reduction coordinates in 2D
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_cell_params
#' @inheritParams plot_cow_params
#' @inheritParams plot_image_params
#' @inheritParams plot_spatenr_params
#' @inheritParams plot_dimred_params
#' @inheritParams plot_nn_net_params
#' @inheritParams plot_params
#' @param plot_alignment direction to align plot
#' @param cell_annotation_values numeric cell annotation columns
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param spat_point_shape shape of points (border, no_border or voronoi)
#' @param spat_point_size size of spatial points
#' @param spat_point_alpha transparency of spatial points
#' @param spat_point_border_col border color of spatial points
#' @param spat_point_border_stroke border stroke of spatial points
#' @param dim_show_cluster_center show the center of each cluster
#' @param dim_show_center_label provide a label for each cluster
#' @param dim_center_point_size size of the center point
#' @param dim_center_point_border_col border color of center point
#' @param dim_center_point_border_stroke stroke size of center point
#' @param dim_label_size size of the center label
#' @param dim_label_fontface font of the center label
#' @param spat_show_cluster_center show the center of each cluster
#' @param spat_show_center_label provide a label for each cluster
#' @param spat_center_point_size size of the spatial center points
#' @param spat_center_point_border_col border color of the spatial center points
#' @param spat_center_point_border_stroke stroke size of the spatial center
#' points
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param dim_edge_alpha column to use for alpha of the edges
#' @param spat_show_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param spat_network_alpha alpha of spatial network
#' @param spat_show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spat_grid_color color of spatial grid
#' @param dim_other_point_size size of not selected dim cells
#' @param spat_other_point_size size of not selected spat cells
#' @param spat_other_cells_alpha alpha of not selected spat cells
#' @param coord_fix_ratio ratio for coordinates
#' @param dim_background_color background color of points in dim. reduction
#' space
#' @param spat_background_color background color of spatial points
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @details Description of parameters.
#' @family spatial and dimension reduction cell annotation visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatDimCellPlot2D(g, cell_annotation_values = "leiden_clus")
#'
#' @export
spatDimCellPlot2D <- function(
        gobject,
        feat_type = NULL,
        spat_unit = NULL,
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        plot_alignment = c("vertical", "horizontal"),
        spat_enr_names = NULL,
        cell_annotation_values = NULL,
        dim_reduction_to_use = "umap",
        dim_reduction_name = "umap",
        dim1_to_use = 1,
        dim2_to_use = 2,
        sdimx = "sdimx",
        sdimy = "sdimy",
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = c("divergent", "sequential"),
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        dim_point_shape = c("border", "no_border"),
        dim_point_size = 1,
        dim_point_alpha = 1,
        dim_point_border_col = "black",
        dim_point_border_stroke = 0.1,
        spat_point_shape = c("border", "no_border", "voronoi"),
        spat_point_size = 1,
        spat_point_alpha = 1,
        spat_point_border_col = "black",
        spat_point_border_stroke = 0.1,
        dim_show_cluster_center = FALSE,
        dim_show_center_label = TRUE,
        dim_center_point_size = 4,
        dim_center_point_border_col = "black",
        dim_center_point_border_stroke = 0.1,
        dim_label_size = 4,
        dim_label_fontface = "bold",
        spat_show_cluster_center = FALSE,
        spat_show_center_label = FALSE,
        spat_center_point_size = 4,
        spat_center_point_border_col = "black",
        spat_center_point_border_stroke = 0.1,
        spat_label_size = 4,
        spat_label_fontface = "bold",
        show_NN_network = FALSE,
        nn_network_to_use = "sNN",
        nn_network_name = "sNN.pca",
        dim_edge_alpha = 0.5,
        spat_show_network = FALSE,
        spatial_network_name = "Delaunay_network",
        spat_network_color = "red",
        spat_network_alpha = 0.5,
        spat_show_grid = FALSE,
        spatial_grid_name = "spatial_grid",
        spat_grid_color = "green",
        show_other_cells = TRUE,
        other_cell_color = "grey",
        dim_other_point_size = 0.5,
        spat_other_point_size = 0.5,
        spat_other_cells_alpha = 0.5,
        show_legend = TRUE,
        legend_text = 8,
        legend_symbol_size = 1,
        dim_background_color = "white",
        spat_background_color = "white",
        vor_border_color = "white",
        vor_max_radius = 200,
        vor_alpha = 1,
        axis_text = 8,
        axis_title = 8,
        coord_fix_ratio = 1,
        cow_n_col = NULL,
        cow_rel_h = 1,
        cow_rel_w = 1,
        cow_align = "h",
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatDimCellPlot2D",
        view = NULL,
        space = NULL) {
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_locs",
            "spatial_enrichment", "expression", "dimension_reduction",
            "images"))
    plot_alignment <- match.arg(plot_alignment,
        choices = c("vertical", "horizontal")
    )

    # dimension reduction plot
    dmpl <- dimCellPlot2D(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        dim_reduction_to_use = dim_reduction_to_use,
        dim_reduction_name = dim_reduction_name,
        dim1_to_use = dim1_to_use,
        dim2_to_use = dim2_to_use,
        spat_enr_names = spat_enr_names,
        cell_annotation_values = cell_annotation_values,
        cell_color_gradient = cell_color_gradient,
        gradient_midpoint = gradient_midpoint,
        gradient_style = gradient_style,
        gradient_limits = gradient_limits,
        select_cell_groups = select_cell_groups,
        select_cells = select_cells,
        point_shape = dim_point_shape,
        point_size = dim_point_size,
        point_alpha = dim_point_alpha,
        point_border_col = dim_point_border_col,
        point_border_stroke = dim_point_border_stroke,
        show_cluster_center = dim_show_cluster_center,
        show_center_label = dim_show_center_label,
        center_point_size = dim_center_point_size,
        center_point_border_col = dim_center_point_border_col,
        center_point_border_stroke = dim_center_point_border_stroke,
        label_size = dim_label_size,
        label_fontface = dim_label_fontface,
        show_NN_network = show_NN_network,
        nn_network_to_use = nn_network_to_use,
        network_name = nn_network_name,
        edge_alpha = dim_edge_alpha,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        other_point_size = dim_other_point_size,
        show_legend = show_legend,
        legend_text = legend_text,
        legend_symbol_size = legend_symbol_size,
        background_color = dim_background_color,
        axis_text = axis_text,
        axis_title = axis_title,
        cow_n_col = cow_n_col,
        cow_rel_h = cow_rel_h,
        cow_rel_w = cow_rel_w,
        cow_align = cow_align,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE
    )

    # spatial plot
    spl <- spatCellPlot2D(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        show_image = show_image,
        gimage = gimage,
        image_name = image_name,
        largeImage_name = largeImage_name,
        sdimx = sdimx,
        sdimy = sdimy,
        spat_enr_names = spat_enr_names,
        cell_annotation_values = cell_annotation_values,
        cell_color_gradient = cell_color_gradient,
        gradient_midpoint = gradient_midpoint,
        gradient_style = gradient_style,
        gradient_limits = gradient_limits,
        select_cell_groups = select_cell_groups,
        select_cells = select_cells,
        point_shape = spat_point_shape,
        point_size = spat_point_size,
        point_alpha = spat_point_alpha,
        point_border_col = spat_point_border_col,
        point_border_stroke = spat_point_border_stroke,
        show_cluster_center = spat_show_cluster_center,
        show_center_label = spat_show_center_label,
        center_point_size = spat_center_point_size,
        center_point_border_col = spat_center_point_border_col,
        center_point_border_stroke = spat_center_point_border_stroke,
        label_size = spat_label_size,
        label_fontface = spat_label_fontface,
        show_network = spat_show_network,
        spatial_network_name = spatial_network_name,
        network_color = spat_network_color,
        network_alpha = spat_network_alpha,
        show_grid = spat_show_grid,
        spatial_grid_name = spatial_grid_name,
        grid_color = spat_grid_color,
        show_other_cells = show_other_cells,
        other_cell_color = other_cell_color,
        other_point_size = spat_other_point_size,
        other_cells_alpha = spat_other_cells_alpha,
        coord_fix_ratio = coord_fix_ratio,
        show_legend = show_legend,
        legend_text = legend_text,
        legend_symbol_size = legend_symbol_size,
        background_color = spat_background_color,
        vor_border_color = vor_border_color,
        vor_max_radius = vor_max_radius,
        vor_alpha = vor_alpha,
        axis_text = axis_text,
        axis_title = axis_title,
        cow_n_col = cow_n_col,
        cow_rel_h = cow_rel_h,
        cow_rel_w = cow_rel_w,
        cow_align = cow_align,
        show_plot = FALSE,
        return_plot = TRUE,
        save_plot = FALSE
    )


    if (plot_alignment == "vertical") {
        ncol <- 1
        nrow <- 2
        combo_plot <- cowplot::plot_grid(
            dmpl, spl,
            ncol = ncol, nrow = nrow,
            rel_heights = c(1), rel_widths = c(1),
            align = "v"
        )
    } else {
        ncol <- 2
        nrow <- 1
        combo_plot <- cowplot::plot_grid(
            dmpl, spl,
            ncol = ncol, nrow = nrow,
            rel_heights = c(1), rel_widths = c(1),
            align = "h"
        )
    }

    return(plot_output_handler(
        gobject = gobject,
        plot_object = combo_plot,
        save_plot = save_plot,
        show_plot = show_plot,
        return_plot = return_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}




#' @title spatDimCellPlot
#' @name spatDimCellPlot
#' @description Visualize numerical features of cells according to spatial
#' AND dimension reduction coordinates in 2D
#' @inheritDotParams spatDimCellPlot2D
#' @details Description of parameters.
#' @family spatial and dimension reduction cell annotation visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' spatDimCellPlot(g, cell_annotation_values = "leiden_clus")
#'
#' @export
spatDimCellPlot <- function(...) {
    spatDimCellPlot2D(...)
}
