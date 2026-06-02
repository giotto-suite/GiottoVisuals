# subcellular ####



#' @title spatInSituPlotPoints
#' @name spatInSituPlotPoints
#' @description Function to plot multiple features for multiple modalities
#' at the spatial in situ level
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_feat_params
#' @inheritParams plot_poly_params
#' @inheritParams plot_image_params
#' @inheritParams plot_spatenr_params
#' @param largeImage_name deprecated
#' @param spat_loc_name name of spatial locations
#' @param feats named list of features to plot
#' @param feat_type feature types of the feats
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param xlim limits of x-scale (min/max vector)
#' @param ylim limits of y-scale (min/max vector)
#' @param remove_background_polygon logical (default = `TRUE`). `xlim` and
#' `ylim` which work through (`crop()`) may
#' sometimes produce extent-filling polygons when the original geometry is
#' problematic or invalid. Set `TRUE` to remove these, based on whether a
#' polygon fills up most of the x and y range.
#' @param point_size size of the points
#' @param stroke stroke to apply to feature points
#' @param expand_counts expand feature coordinate counts (see details)
#' @param count_info_column column name with count information
#' (if expand_counts = TRUE)
#' @param jitter numeric. Maximum x,y jitter provided as c(x, y) or a single
#' number which will be recycled to length 2.
#' @param axis_text axis text size
#' @param axis_title title text size
#' @param legend_text legend text size
#' @param coord_fix_ratio fix ratio of coordinates
#' @param background_color background color
#' @param show_legend show legend
#' @param plot_method method to plot points
#' @param plot_last which layer to show on top of plot,
#' polygons (default) or points.
#' @param theme_param list of additional params passed to `ggplot2::theme()`
#' @param verbose be verbose
#' @param view,space optional [GiottoClass::giottoView-class] /
#' [GiottoClass::giottoSpace-class] or the name of one slotted on
#' `gobject`. When supplied, `gobject` is pre-narrowed via
#' [GiottoClass::materialize()] before any sub-layer fetches its data.
#' One resolver pass is shared across the polygon / point / image
#' helpers. Backend-agnostic: in-mem and on-disk gobjects take the
#' same code path via the registered view coordinator.
#' @param ... additional params to pass to [spatValues()]
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("vizgen")
#'
#' # plot the number detected features in the `giotto` object.
#' spatInSituPlotPoints(g,
#'     polygon_feat_type = "aggregate",
#'     polygon_fill = "nr_feats",
#'     polygon_fill_as_factor = FALSE,
#'     polygon_fill_gradient_style = "sequential",
#'     polygon_line_size = 0.1,
#'     polygon_alpha = 1
#' )
#'
#' # plot the same as above, but with the first 4 rna features plotted as
#' # detection points. Also add in the background and change the polygon
#' # alpha and border color
#' spatInSituPlotPoints(g,
#'     polygon_feat_type = "aggregate",
#'     polygon_fill = "nr_feats",
#'     polygon_fill_as_factor = FALSE,
#'     polygon_fill_gradient_style = "sequential",
#'     polygon_line_size = 0.1,
#'     polygon_alpha = 0.4,
#'     polygon_color = "magenta",
#'     feats = list("rna" = GiottoClass::featIDs(g)[1:4]),
#'     point_size = 0.8,
#'     plot_last = "points",
#'     show_image = TRUE,
#'     image_name = "dapi_z0"
#' )
#'
#' # plot with spatial enrichment information
#' spatInSituPlotPoints(g,
#'     polygon_feat_type = "aggregate",
#'     spat_enr_name = "cluster_metagene",
#'     polygon_fill = "1",
#'     polygon_fill_as_factor = FALSE,
#'     polygon_fill_gradient_style = "sequential",
#'     polygon_alpha = 1
#' )
#'
#' # plot with xlim,ylim and feature expression
#' spatInSituPlotPoints(g,
#'     polygon_feat_type = "aggregate",
#'     xlim = c(6500, 6920),
#'     ylim = c(-5100, -4820),
#'     polygon_alpha = 1,
#'     polygon_fill = "Slc17a7",
#'     polygon_fill_gradient_midpoint = 0,
#'     expression_values = "scaled",
#'     show_image = TRUE
#' )
#'
#' @family In Situ visualizations
#' @export
spatInSituPlotPoints <- function(
        gobject,
        show_image = FALSE,
        gimage = NULL,
        image_name = NULL,
        largeImage_name = NULL,
        spat_unit = NULL,
        spat_loc_name = NULL,
        feats = NULL,
        feat_type = "rna",
        feats_color_code = NULL,
        feat_shape_code = NULL,
        sdimx = "x",
        sdimy = "y",
        xlim = NULL,
        ylim = NULL,
        remove_background_polygon = TRUE,
        spat_enr_name = NULL,
        point_size = 1.5,
        stroke = 0.5,
        expand_counts = FALSE,
        count_info_column = "count",
        jitter = c(0, 0),
        show_polygon = TRUE,
        use_overlap = TRUE,
        polygon_feat_type = "cell",
        polygon_color = "grey",
        polygon_bg_color = "black",
        polygon_fill = NULL,
        polygon_fill_gradient = NULL,
        polygon_fill_gradient_midpoint = NULL,
        polygon_fill_gradient_style = c("divergent", "sequential"),
        polygon_fill_as_factor = NULL,
        polygon_fill_code = NULL,
        polygon_alpha = NULL,
        polygon_line_size = 0.4,
        axis_text = 8,
        axis_title = 8,
        legend_text = 6,
        coord_fix_ratio = 1,
        background_color = "black",
        show_legend = TRUE,
        show_axes = NULL,
        plot_method = c("ggplot", "scattermore", "scattermost"),
        plot_last = c("polygons", "points"),
        theme_param = list(),
        show_plot = NULL,
        return_plot = NULL,
        save_plot = NULL,
        save_param = list(),
        default_save_name = "spatInSituPlotPoints",
        verbose = TRUE,
        view = NULL,
        space = NULL,
        ...) {
    if (inherits(gobject, "giottoMulti")) {
        return(.gg_multi_dispatch_spatial(
            plot_fn = spatInSituPlotPoints,
            named = mget(setdiff(names(formals()), "...")),
            dots = list(...),
            gobject = gobject, view = view, space = space
        ))
    }
    .gg_assert_giotto_single(gobject)
    # currently not intended for more than one thing to plot
    checkmate::assert_character(polygon_fill, null.ok = TRUE, len = 1L)
    handle_errors({
    # Pre-narrow once for the slots this plot reads (polygons, points,
    # locations, enrichment, expression-for-fill, metadata, images).
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_info", "spatial_locs",
            "spatial_enrichment", "feat_info", "feat_metadata",
            "expression", "images"))
    # set polygon_feat_type
    avail_poly_names <- list_spatial_info_names(gobject = gobject)
    if (polygon_feat_type == "cell" &&
        !"cell" %in% avail_poly_names) {
        polygon_feat_type <- spat_unit
        if (verbose) {
            wrap_msg(
                "[polygon_feat_type] 'cell' not discovered in polygon names.
            Defaulting to spat_unit."
            )
        }
    }

    send_warn <- getOption("giotto.warn_sispp_feats", TRUE)
    if (is.null(feats) && send_warn) {
        warning(wrap_txt(
            "You need to select features (feats) and modify feature
            types (feat_type) if you want to show individual features
            (e.g. transcripts)
            This warning is shown once per session"
        ))
        options("giotto.warn_sispp_feats" = FALSE)
    }

    # deprecation message
    if (!is.null(largeImage_name)) {
        deprecate_warn(
            when = "0.2.0",
            what = "spatInSituPlotPoints(largeImage_name)",
            details = c(
                "Use `image_name` argument instead for all images to plot."
            )
        )
        image_name <- c(image_name, largeImage_name)
    }

    # check valid input
    plot_last <- match.arg(plot_last, choices = c("polygons", "points"))

    ## giotto image ##
    if (isTRUE(show_image)) {
        # get 1 or more images
        gimage <- getGiottoImage(
            gobject = gobject,
            name = image_name
        )
    }

    if (isTRUE(show_image)) {
        polygon_alpha <- polygon_alpha %null% 0.5
    } else {
        polygon_alpha <- polygon_alpha %null% 1
    }

    show_axes <- show_axes %null% TRUE

    # what to plot...
    plot_order <- character(0L)
    plot_order <- c(plot_order, "image"[isTRUE(show_image) && !is.null(gimage)])
    plot_order <- switch(plot_last,
        "polygons" = c(plot_order,
            "points"[!is.null(feats)], "polygons"[show_polygon]
        ),
        "points" = c(plot_order,
            "polygons"[show_polygon], "points"[!is.null(feats)]
        )
    )
    if (length(plot_order) == 0L) {
        vmsg(.v = verbose, "No elements to plot")
        return(invisible())
    }


    # plotting funs
    # preload with params to make order switching easy
    img_fun <- function(p) {
        .sissp_image(p,
            gobject = gobject,
            xlim = xlim,
            ylim = ylim,
            spat_unit = spat_unit,
            feat_type = feat_type,
            spat_loc_name = spat_loc_name,
            polygon_feat_type = polygon_feat_type,
            gimage = gimage,
            verbose = verbose
        )
    }
    poly_fun <- function(p) {
        .sissp_polygon(p,
            gobject = gobject,
            polygon_feat_type = polygon_feat_type,
            feat_type = feat_type,
            remove_background_polygon = remove_background_polygon,
            xlim = xlim,
            ylim = ylim,
            polygon_fill = polygon_fill,
            spat_loc_name = spat_loc_name,
            spat_enr_name = spat_enr_name,
            show_legend = show_legend,
            polygon_fill_gradient = polygon_fill_gradient,
            polygon_fill_gradient_midpoint = polygon_fill_gradient_midpoint,
            polygon_fill_gradient_style = polygon_fill_gradient_style,
            polygon_fill_as_factor = polygon_fill_as_factor,
            polygon_fill_code = polygon_fill_code,
            polygon_bg_color = polygon_bg_color,
            polygon_color = polygon_color,
            polygon_alpha = polygon_alpha,
            polygon_line_size = polygon_line_size,
            verbose = verbose,
            ...
        )
    }
    pt_fun <- function(p) {
        .sissp_points(p,
            gobject = gobject,
            feats = feats,
            sdimx = sdimx,
            sdimy = sdimy,
            polygon_feat_type = polygon_feat_type,
            feat_type = feat_type,
            use_overlap = use_overlap,
            feats_color_code = feats_color_code,
            feat_shape_code = feat_shape_code,
            expand_counts = expand_counts,
            count_info_column = count_info_column,
            jitter = jitter,
            point_size = point_size,
            stroke = stroke,
            show_legend = show_legend,
            plot_method = plot_method,
            verbose = verbose
        )
    }

    ## 1. start plot
    plot <- ggplot2::ggplot()

    ## 2. add plot elements
    for (elem in plot_order) {
        elem_fun <- switch(elem,
            "image" = img_fun,
            "polygons" = poly_fun,
            "points" = pt_fun
        )
        plot <- elem_fun(plot)
    }

    ## 3. adjust theme settings
    gg_theme_args <- c(
        theme_param,
        legend_text = legend_text,
        axis_title = axis_title,
        axis_text = axis_text,
        background_color = background_color
    )
    plot <- plot + do.call(.gg_theme, args = gg_theme_args)
    if (!show_axes) {
        plot <- .theme_remove_axes(plot)
    }

    # subset data based on x and y limits
    if (!is.null(xlim)) {
        plot <- plot + ggplot2::xlim(xlim)
    }
    if (!is.null(ylim)) {
        plot <- plot + ggplot2::ylim(ylim)
    }

    # fix coordinates
    plot <- .aspect_ratio(plot, coord_fix_ratio)

    return(plot_output_handler(
        gobject = gobject,
        plot_object = plot,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
    }, prefix = "spatInSituPlotPoints()")
}


.sissp_image <- function(plot, gobject,
    xlim = NULL,
    ylim = NULL,
    spat_unit = NULL,
    feat_type = NULL,
    spat_loc_name = NULL,
    polygon_feat_type = NULL,
    gimage = NULL,
    verbose = NULL
) {
    e <- NULL
    if (!is.null(xlim) && !is.null(ylim)) {
        e <- ext(c(xlim, ylim))
    }

    plot <- plot_spat_image_layer_ggplot(
        gg_obj = plot,
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        spat_loc_name = spat_loc_name,
        polygon_feat_type = polygon_feat_type,
        gimage = gimage,
        ext = e
    )

    vmsg(.v = verbose, "plot image layer done")
    plot
}

# Record `xlim` / `ylim` as a crop step on a scratch view, and return the
# gobject carrying it together with that view's name.
#
# The window used to be handed to `combineCellData()`, which cropped the
# polygon carrier eagerly. That works on a `SpatVector` and not on a
# disk-backed store, where neither `crop()` nor `terra::geom()` has a route --
# a crop STEP does, because it resolves through whichever view coordinator the
# object registers. So the window is expressed the same way a user-supplied
# view is, and takes the same path.
#
# The view is scratch: the gobject carrying it is local to the plotting call
# and never returned, so nothing is left on the user's object. The name is
# chosen to avoid collision with anything already slotted.
#
# `geom = "poly"` rather than the `centroid` default, to match the eager
# behaviour this replaces -- cropping the polygon carrier kept a polygon
# overlapping the window even when its centroid fell outside.
#
# Stacking: a caller-supplied view is copied onto the scratch name first, so
# the window lands last. That ordering wants revisiting when `view` / `space`
# become routine here; it is the reading that matches what the eager crop did.
.sissp_window_view <- function(gobject, spat_unit,
                               xlim = NULL, ylim = NULL, view = NULL) {
    if (is.null(xlim) && is.null(ylim)) {
        return(list(gobject = gobject, view = view))
    }

    # axes the caller left open keep the object's own extent, which is what
    # starting from the carrier's extent used to accomplish
    e <- terra::ext(GiottoClass::ext(gobject,
        spat_unit = spat_unit, prefer = "polygon"))
    if (!is.null(xlim)) e[c(1, 2)] <- xlim
    if (!is.null(ylim)) e[c(3, 4)] <- ylim

    # giottoViews() returns the names themselves, not a named list
    taken <- GiottoClass::giottoViews(gobject)
    nm <- ".roi"
    i <- 1L
    while (nm %in% taken) {
        nm <- sprintf(".roi%d", i)
        i <- i + 1L
    }

    if (!is.null(view)) {
        GiottoClass::giottoView(gobject, nm) <-
            GiottoClass::giottoView(gobject, view)
    }
    gobject <- GiottoClass::crop(gobject, e, geom = "poly", view = nm)
    list(gobject = gobject, view = nm)
}

.sissp_polygon <- function(plot, gobject,
    polygon_feat_type = NULL,
    feat_type = NULL,
    remove_background_polygon = TRUE,
    xlim = NULL,
    ylim = NULL,
    polygon_fill = NULL,
    spat_loc_name = NULL,
    spat_enr_name = NULL,
    show_legend = TRUE,
    polygon_fill_gradient = NULL,
    polygon_fill_gradient_midpoint = NULL,
    polygon_fill_gradient_style = c("divergent", "sequential"),
    polygon_fill_as_factor = NULL,
    polygon_fill_code = NULL,
    polygon_bg_color = "black",
    polygon_color = "grey",
    polygon_alpha = NULL,
    polygon_line_size = 0.4,
    verbose = NULL,
    ...
) {
    # Set feat_type and spat_unit
    polygon_feat_type <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = polygon_feat_type
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = polygon_feat_type,
        feat_type = feat_type
    )

    # the window travels as a view, not as a carrier crop -- see
    # .sissp_window_view() for why
    .win <- .sissp_window_view(gobject, spat_unit = polygon_feat_type,
        xlim = xlim, ylim = ylim)

    polygon_combo <- combineCellData(
        gobject = .win$gobject,
        spat_loc_name = spat_loc_name,
        feat_type = feat_type,
        include_poly_info = TRUE,
        poly_info = polygon_feat_type,
        remove_background_polygon = remove_background_polygon,
        view = .win$view
    )

    polygon_dt <- data.table::rbindlist(polygon_combo, fill = TRUE)

    if (!is.null(polygon_fill) &&
        !polygon_fill %in% names(polygon_dt)) {
        data <- spatValues(gobject,
            feats = polygon_fill,
            spat_unit = polygon_feat_type,
            feat_type = feat_type,
            spat_enr_name = spat_enr_name,
            ...
        )
        polygon_dt <- merge(polygon_dt, data,
            by = "cell_ID", all.x = TRUE)
    }

    data.table::setnames(polygon_dt, old = "cell_ID", new = "poly_ID")

    plot <- plot_cell_polygon_layer(
        ggobject = plot,
        instrs = instructions(gobject),
        polygon_dt = polygon_dt,
        polygon_grouping = "poly_ID",
        fill = polygon_fill,
        poly_fill_gradient = polygon_fill_gradient,
        fill_gradient_midpoint = polygon_fill_gradient_midpoint,
        fill_gradient_style = polygon_fill_gradient_style,
        fill_as_factor = polygon_fill_as_factor,
        fill_code = polygon_fill_code,
        bg_color = polygon_bg_color,
        color = polygon_color,
        alpha = polygon_alpha,
        linewidth = polygon_line_size,
        show_legend = show_legend
    )

    vmsg(.v = verbose, "plot polygon layer done")
    plot
}

.sissp_points <- function(plot, gobject,
    feats = NULL,
    sdimx = "x",
    sdimy = "y",
    polygon_feat_type = NULL,
    feat_type = NULL,
    use_overlap = TRUE,
    feats_color_code = NULL,
    feat_shape_code = NULL,
    expand_counts = FALSE,
    count_info_column = "count",
    jitter = c(0, 0),
    point_size = 1.5,
    stroke = 0.5,
    show_legend = TRUE,
    plot_method = c("ggplot", "scattermore", "scattermost"),
    verbose = NULL
) {
    # use_overlap = TRUE will use the overlap results
    # use_overlap = FALSE will use the raw/full tx coordinate results
    if (isTRUE(use_overlap)) {
        # TODO: check if overlap exists, if not print warning message
        # and default to non-overlap results
        spatial_feat_info <- combineFeatureOverlapData(
            gobject = gobject,
            feat_type = feat_type,
            sel_feats = feats,
            poly_info = polygon_feat_type
        )
    } else {
        spatial_feat_info <- combineFeatureData(
            gobject = gobject,
            spat_unit = polygon_feat_type,
            feat_type = feat_type,
            sel_feats = feats
        )
    }

    spatial_feat_info <- data.table::rbindlist(spatial_feat_info, fill = TRUE)

    plot <- plot_feature_points_layer(
        ggobject = plot,
        instrs = instructions(gobject),
        ext = ext(gobject, prefer = "points"),
        spatial_feat_info = spatial_feat_info,
        feats = feats,
        feats_color_code = feats_color_code,
        feat_shape_code = feat_shape_code,
        expand_counts = expand_counts,
        count_info_column = count_info_column,
        jitter = jitter,
        sdimx = sdimx,
        sdimy = sdimy,
        color = "feat_ID",
        shape = "feat",
        point_size = point_size,
        stroke = stroke,
        show_legend = show_legend,
        plot_method = plot_method
    )

    vmsg(.v = verbose, "plot feature points layer done")
    plot
}












# hexbin ####

#' @title Spatial in-situ hexbin plot - single
#' @name .spatInSituPlotHex_single
#' @description function to plot hexbins at the spatial in situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
.spatInSituPlotHex_single <- function(gobject,
    feat = NULL,
    feat_type = "rna",
    sdimx = "x",
    sdimy = "y",
    binwidth = NULL,
    min_axis_bins = NULL,
    alpha = 0.5,
    show_polygon = TRUE,
    polygon_feat_type = "cell",
    polygon_color = "black",
    polygon_fill = NULL,
    polygon_fill_as_factor = NULL,
    polygon_alpha = 0.5,
    polygon_size = 0.5,
    coord_fix_ratio = NULL,
    axis_text = 8,
    axis_title = 8,
    legend_text = 6,
    background_color = "black") {
    if (is.null(feat)) {
        stop("You need to select a feature (feat) and modify feature
            types (feat_type) if needed \n")
    }

    plot <- ggplot2::ggplot()

    ## polygon layer ##
    if (show_polygon == TRUE) {
        if (is.null(polygon_feat_type)) {
            polygon_feat_type <- gobject@expression_feat[[1]]
        }

        polygon_dt <- getPolygonInfo(
            gobject = gobject,
            name = polygon_feat_type
        ) %>%
            data.table::as.data.table(geom = "XY")

        plot <- plot_cell_polygon_layer(
            ggobject = plot,
            instrs = instructions(gobject),
            polygon_dt,
            polygon_grouping = "poly_ID",
            fill = polygon_fill,
            fill_as_factor = polygon_fill_as_factor,
            color = polygon_color,
            alpha = polygon_alpha,
            linewidth = polygon_size
        )
    }



    ## hexbin layer ##

    form_feat <- list(feat_type = c(feat))
    spatial_feat_info <- combineFeatureOverlapData(
        gobject = gobject,
        feat_type = feat_type,
        sel_feats = form_feat,
        poly_info = polygon_feat_type
    )

    spatial_feat_info <- do.call("rbind", spatial_feat_info)

    plot <- plot_feature_hexbin_layer(
        ggobject = plot,
        instrs = instructions(gobject),
        spatial_feat_info = spatial_feat_info,
        sel_feat = feat,
        sdimx = sdimx,
        sdimy = sdimy,
        binwidth = binwidth,
        min_axis_bins = min_axis_bins
    )


    ## adjust theme settings
    plot <- plot + ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_text),
        axis.title = element_text(size = axis_title),
        axis.text = element_text(size = axis_text),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = background_color)
    )

    # fix coord ratio
    plot <- .aspect_ratio(plot, coord_fix_ratio)

    return(plot)
}



#' @title spatInSituPlotHex
#' @name spatInSituPlotHex
#' @description Function to plot hexbins for features for multiple
#' modalities at the spatial in situ level
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_poly_params
#' @inheritParams plot_cow_params
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param binwidth numeric vector for x and y width of bins
#' (default is minor axis range/10, where the 10 is from \code{min_axis_bins})
#' @param min_axis_bins number of bins to create per range defined by minor axis
#' (default value is 10)
#' @param alpha alpha of hexbin plot
#' @param polygon_size deprecated
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param axis_text axis text size
#' @param axis_title title text size
#' @param legend_text legend text size
#' @param background_color background color
#' @details TODO
#' @family In Situ visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("vizgen")
#' spatInSituPlotHex(g,
#'     feats = c("Mlc1", "Gprc5b", "Gfap"),
#'     polygon_feat_type = "z0"
#' )
#' @export
spatInSituPlotHex <- function(gobject,
    feats = NULL,
    feat_type = "rna",
    sdimx = "x",
    sdimy = "y",
    binwidth = NULL,
    min_axis_bins = 10,
    alpha = 0.5,
    show_polygon = TRUE,
    polygon_feat_type = "cell",
    polygon_color = "black",
    polygon_fill = NULL,
    polygon_fill_as_factor = NULL,
    polygon_alpha = 0.5,
    polygon_size = deprecated(),
    polygon_line_size = 0.5,
    coord_fix_ratio = 1,
    axis_text = 8,
    axis_title = 8,
    legend_text = 6,
    background_color = "white",
    cow_n_col = NULL,
    cow_rel_h = 1,
    cow_rel_w = 1,
    cow_align = "h",
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "spatInSituPlotHex",
    view = NULL,
    space = NULL) {
    if (inherits(gobject, "giottoMulti")) {
        return(.gg_multi_dispatch_spatial(
            plot_fn = spatInSituPlotHex,
            named = mget(names(formals())),
            gobject = gobject, view = view, space = space
        ))
    }
    .gg_assert_giotto_single(gobject)
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "spatial_info", "spatial_locs",
            "feat_info", "feat_metadata", "images"))
    # deprecate
    if (GiottoUtils::is_present(polygon_size)) {
        deprecate_warn(
            "0.0.0.9000",
            "GiottoVisuals::spatInSituPlotHex(polygon_size = )",
            "GiottoVisuals::spatInSituPlotHex(polygon_line_size = )"
        )
        polygon_line_size <- polygon_size
    }

    package_check("hexbin")

    if (is.null(feats)) {
        stop("You need to select features (feats) and modify
            feature types (feat_type) if needed \n")
    }

    # print, return and save parameters
    show_plot <- ifelse(is.null(show_plot),
        readGiottoInstructions(gobject,
            param = "show_plot"
        ), show_plot
    )
    save_plot <- ifelse(is.null(save_plot),
        readGiottoInstructions(gobject,
            param = "save_plot"
        ), save_plot
    )
    return_plot <- ifelse(is.null(return_plot),
        readGiottoInstructions(gobject,
            param = "return_plot"
        ),
        return_plot
    )

    ## plotting ##
    savelist <- list()

    for (sel_feat in feats) {
        pl <- .spatInSituPlotHex_single(
            gobject = gobject,
            feat = sel_feat,
            feat_type = feat_type,
            sdimx = sdimx,
            sdimy = sdimy,
            binwidth = binwidth,
            min_axis_bins = min_axis_bins,
            alpha = alpha,
            show_polygon = show_polygon,
            polygon_feat_type = polygon_feat_type,
            polygon_color = polygon_color,
            polygon_fill = polygon_fill,
            polygon_fill_as_factor = polygon_fill_as_factor,
            polygon_alpha = polygon_alpha,
            polygon_size = polygon_line_size,
            coord_fix_ratio = coord_fix_ratio,
            axis_text = axis_text,
            axis_title = axis_title,
            legend_text = legend_text,
            background_color = background_color
        )

        savelist[[sel_feat]] <- pl
    }

    if (length(savelist) == 1) {
        combo_plot <- savelist[[1]]
    } else {
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
    }


    ## print plot
    if (show_plot == TRUE) {
        print(combo_plot)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call("all_plots_save_function", c(
            list(
                gobject = gobject,
                plot_object = combo_plot,
                default_save_name = default_save_name
            ),
            save_param
        ))
    }

    ## return plot
    if (return_plot == TRUE) {
        return(combo_plot)
    }
}








# density ####

#' @title Spatial in-situ density plot - single
#' @name .spatInSituPlotDensity_single
#' @description low level function to plot density plots at the spatial in
#' situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
.spatInSituPlotDensity_single <- function(gobject,
    data,
    feat = NULL,
    sdimx = "x",
    sdimy = "y",
    alpha = 0.95,
    show_polygon = TRUE,
    polygon_feat_type = "cell",
    polygon_color = "black",
    polygon_fill = NULL,
    polygon_fill_as_factor = NULL,
    polygon_alpha = 0.5,
    polygon_size = 0.5,
    coord_fix_ratio = NULL,
    axis_text = 8,
    axis_title = 8,
    legend_text = 6,
    background_color = "black") {

    feat_ID <- NULL # NSE vars

    if (is.null(feat)) {
        stop("You need to select a feature (feat) and modify feature types
            (feat_type) if needed \n")
    }

    plot <- ggplot2::ggplot()

    feat_data <- data[feat_ID == feat,]

    plot <- plot_feature_raster_density_layer(
        ggobject = plot,
        instrs = instructions(gobject),
        spatial_feat_info = feat_data,
        sel_feat = feat,
        sdimx = sdimx,
        sdimy = sdimy,
        alpha = alpha
    )

    ## polygon layer ##
    if (show_polygon == TRUE) {
        if (is.null(polygon_feat_type)) {
            polygon_feat_type <- gobject@expression_feat[[1]]
        }


        polygon_dt <- getPolygonInfo(
            gobject = gobject,
            name = polygon_feat_type
        ) %>%
            data.table::as.data.table(geom = "XY")

        plot <- plot_cell_polygon_layer(
            ggobject = plot,
            instrs = instructions(gobject),
            polygon_dt,
            polygon_grouping = "poly_ID",
            fill = polygon_fill,
            fill_as_factor = polygon_fill_as_factor,
            color = polygon_color,
            alpha = polygon_alpha,
            linewidth = polygon_size
        )
    }


    ## adjust theme settings
    plot <- plot + ggplot2::theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = legend_text),
        axis.title = element_text(size = axis_title),
        axis.text = element_text(size = axis_text),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = background_color)
    )

    # fix coord ratio
    plot <- .aspect_ratio(plot, coord_fix_ratio)

    return(plot)
}



#' @title spatInSituPlotDensity
#' @name spatInSituPlotDensity
#' @description Function for density plots for features for multiple
#' modalities at the spatial in situ level
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @inheritParams plot_poly_params
#' @inheritParams plot_cow_params
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param alpha alpha of density plot
#' @param use_overlap whether features should be restricted to those overlapped
#' by the `polygon_feat_type` spatial unit.
#' @param polygon_size deprecated
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param axis_text axis text size
#' @param axis_title title text size
#' @param legend_text legend text size
#' @param background_color background color
#' @details TODO
#' @family In Situ visualizations
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("vizgen")
#' spatInSituPlotDensity(g,
#'     feats = c("Mlc1", "Gprc5b", "Gfap"),
#'     polygon_feat_type = "z0"
#' )
#'
#' @export
spatInSituPlotDensity <- function(gobject,
    feats = NULL,
    feat_type = "rna",
    sdimx = "x",
    sdimy = "y",
    alpha = 0.95,
    show_polygon = TRUE,
    use_overlap = FALSE,
    polygon_feat_type = "cell",
    polygon_color = "black",
    polygon_fill = NULL,
    polygon_fill_as_factor = NULL,
    polygon_alpha = 0.5,
    polygon_size = deprecated(),
    polygon_line_size = 0.5,
    coord_fix_ratio = 1,
    axis_text = 8,
    axis_title = 8,
    legend_text = 6,
    background_color = "black",
    cow_n_col = NULL,
    cow_rel_h = 1,
    cow_rel_w = 1,
    cow_align = "h",
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "spatInSituPlotDensity",
    view = NULL,
    space = NULL) {
    if (inherits(gobject, "giottoMulti")) {
        return(.gg_multi_dispatch_spatial(
            plot_fn = spatInSituPlotDensity,
            named = mget(names(formals())),
            gobject = gobject, view = view, space = space
        ))
    }
    .gg_assert_giotto_single(gobject)
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("spatial_info", "spatial_locs",
            "feat_info", "feat_metadata", "images"))
    # deprecate
    if (GiottoUtils::is_present(polygon_size)) {
        deprecate_warn(
            "0.0.0.9000",
            "GiottoVisuals::spatInSituPlotDensity(polygon_size = )",
            "GiottoVisuals::spatInSituPlotDensity(polygon_line_size = )"
        )
        polygon_line_size <- polygon_size
    }


    if (is.null(feats)) {
        stop("You need to select features (feat) and modify feature types
            (feat_type) if needed \n")
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

    if (isTRUE(use_overlap)) {
        feats2get <- list(feat_type = feats)
        data <- combineFeatureOverlapData(gobject,
            feat_type = feat_type,
            sel_feats = feats2get,
            poly_info = polygon_feat_type
        )
        data <- do.call(rbind, data)
    } else {
        data <- getFeatureInfo(gobject,
            feat_type = feat_type,
            return_giottoPoints = TRUE
        )
        data <- data[feats]
        data <- data.table::as.data.table(data, geom = "XY")
    }

    data <- data[, c("feat_ID", "x", "y")]

    ## plotting ##
    savelist <- list()

    for (sel_feat in feats) {
        pl <- .spatInSituPlotDensity_single(
            gobject = gobject,
            data = data,
            feat = sel_feat,
            sdimx = sdimx,
            sdimy = sdimy,
            alpha = alpha,
            show_polygon = show_polygon,
            polygon_feat_type = polygon_feat_type,
            polygon_color = polygon_color,
            polygon_fill = polygon_fill,
            polygon_fill_as_factor = polygon_fill_as_factor,
            polygon_alpha = polygon_alpha,
            polygon_size = polygon_line_size,
            coord_fix_ratio = coord_fix_ratio,
            axis_text = axis_text,
            axis_title = axis_title,
            legend_text = legend_text,
            background_color = background_color
        )

        savelist[[sel_feat]] <- pl
    }

    if (length(savelist) == 1) {
        combo_plot <- savelist[[1]]
    } else {
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
    }



    ## print plot
    if (show_plot == TRUE) {
        print(combo_plot)
    }

    ## save plot
    if (save_plot == TRUE) {
        do.call("all_plots_save_function", c(
            list(
                gobject = gobject,
                plot_object = combo_plot,
                default_save_name = default_save_name
            ),
            save_param
        ))
    }

    ## return plot
    if (return_plot == TRUE) {
        return(combo_plot)
    }
}






# helpers ####







## other ####

#' @title expand_feature_info
#' @name expand_feature_info
#' @description low level function to expand feature coordinates by adding
#' jitter to coordinates
#' @param spatial_feat_info a data.table of spatial feature information
#' @param expand_counts logical. whether points should be expanded based on
#' count
#' @param jitter numeric vector. amount of jitter to add
#' @param verbose be verbose
#' @returns data.table
#' @examples
#' g <- GiottoData::loadGiottoMini("vizgen")
#' x <- GiottoClass::getFeatureInfo(g)
#' x <- data.table::as.data.table(x)
#'
#' @keywords internal
#' @export
expand_feature_info <- function(spatial_feat_info,
    expand_counts = FALSE,
    count_info_column = "count",
    jitter = c(0, 0),
    verbose = TRUE) {
    # data.table variables
    feat_ID <- x <- y <- feat <- spat_unit <- NULL

    # 1. expand feature locations with multiple counts
    # (e.g. in seq-Scope or Stereo-seq)
    if (isTRUE(expand_counts)) {
        if (!count_info_column %in% colnames(spatial_feat_info)) {
            stop("count_info_column ", count_info_column, " does not exist")
        }

        if (isTRUE(verbose)) {
            wrap_msg("Start expanding feature information based on count
                    column")
        }


        extra_feats <- spatial_feat_info[get(count_info_column) > 1]
        extra_feats <- extra_feats[, rep(
            get(count_info_column),
            get(count_info_column)
        ),
        by = .(feat_ID, x, y, feat, spat_unit)
        ]
        spatial_feat_info <- rbind(
            extra_feats[
                ,
                .(feat_ID, x, y, feat, spat_unit)
            ],
            spatial_feat_info[
                get(count_info_column) == 1,
                .(feat_ID, x, y, feat, spat_unit)
            ]
        )
    }

    # 2. add jitter to x and y coordinates

    if (!identical(c(0, 0), jitter)) {
        if (isTRUE(verbose)) {
            wrap_msg("Start adding jitter to x and y based on provided
                    max jitter information")
        }

        # create jitter for x and y coordinates: from 0 to max-x or max-y
        tx_number <- nrow(spatial_feat_info)
        x_jitter <- sample(0:jitter[[1]], size = tx_number, replace = TRUE)
        y_jitter <- sample(0:jitter[[2]], size = tx_number, replace = TRUE)

        spatial_feat_info[, c("x", "y") := list(x + x_jitter, y + y_jitter)]
    }

    return(spatial_feat_info)
}
