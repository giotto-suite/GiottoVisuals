#' @include aux_save.R
NULL


# view / space pre-narrow ####

# Thin wrapper around `GiottoClass::materialize()` for plot functions.
#
# Why this helper exists (don't inline it back):
# `GiottoClass::materialize()` dispatches on `signature(giotto,
# giottoView)` / `(giotto, character)` / multi variants — there is no
# method for `view = NULL`. Calling `materialize(g, NULL, NULL)` would
# error with "no applicable method". Plot functions need to accept
# `view = NULL` / `space = NULL` as the default (so callers without a
# view see normal behavior), so each call site would otherwise need an
# inline `if (!is.null(view) || !is.null(space)) ...` guard. This
# helper folds that guard into one place.
#
# `slots` is a character vector of the slot names this plot reads —
# see [GiottoClass::materialize()] for the canonical set. The resolver
# runs ONCE for the listed slots; the predicate / crop is evaluated
# once regardless of how many combine* / getter calls the plot
# internals make.
#
# Backend-agnostic by design: `materialize()` auto-selects
# `dataTableCoordinator` for in-mem gobjects and `parquetCoordinator`
# for `gsource`-backed gobjects via `.default_view_coordinator`.
#
#' @keywords internal
#' @noRd
.gg_materialize <- function(gobject, view, space, slots) {
    if (is.null(view) && is.null(space)) return(gobject)
    GiottoClass::materialize(gobject, view, space = space, slots = slots)
}


# coord fixed ratio ####

.aspect_ratio <- function(pl, coord_fix_ratio = NULL) {
    if (is.null(coord_fix_ratio)) return(pl)
    if (identical(coord_fix_ratio, TRUE)) coord_fix_ratio <- 1
    pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
}

# groupby ####

# split apart a table of data based on a group by column
.groupby <- function(data, group_by, group_by_subset = NULL) {
    possible_groups <- colnames(data)

    ## error if group_by col is not found
    if (!group_by %in% possible_groups) {
        stop("group_by ", group_by, " was not found in pDataDT()")
    }

    # subset unique_groups if needed
    if (!is.null(group_by_subset)) {
        unique_groups <- unique(data[[group_by]])
        not_found <- group_by_subset[!group_by_subset %in% unique_groups]
        if (length(not_found) > 0) {
            message("the following subset was not found: ", not_found)
        }
        unique_groups <- unique_groups[unique_groups %in% group_by_subset]

        # subset data to only group_by_subset
        data <- data[get(group_by) %in% unique_groups]
    }

    datalist <- split(data, data[[group_by]])
}


# clusters ####

#' @title Decide cluster order
#' @name .decide_cluster_order
#' @description creates order for clusters
#' @inheritParams data_access_params
#' @param expression_values expression values to use
#' (e.g. "normalized", "scaled", "custom")
#' @param feats features to use (e.g. genes)
#' @param cluster_column name of column to use for clusters
#' (e.g. "leiden_clus")
#' @param cluster_order method to determine cluster order
#' (e.g. "size", "correlation", "custom")
#' @param cluster_custom_order custom order for clusters
#' @param cor_method method for correlation, default to 'pearson'
#' @param hclust_method method for hierarchical clustering, default to 'ward.D'
#' @return custom
#' @details Calculates order for clusters.
#' @keywords internal
.decide_cluster_order <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    feats,
    cluster_column = NULL,
    cluster_order = c("size", "correlation", "custom"),
    cluster_custom_order = NULL,
    cor_method = "pearson",
    hclust_method = "ward.D") {
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

    # expression data
    values <- match.arg(
        expression_values,
        unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
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
    # subset_values <- expr_values[rownames(expr_values) %in% detected_feats, ]

    # metadata
    cell_metadata <- pDataDT(gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    ## check parameters
    if (is.null(cluster_column)) .gstop("cluster column must be selected")
    if (!cluster_column %in% colnames(cell_metadata)) {
        .gstop(
            "cluster column is not found in",
            str_bracket(spat_unit), str_bracket(feat_type),
            "metadata"
        )
    }

    ## cluster order ##
    cluster_order <- match.arg(
        cluster_order,
        c("size", "correlation", "custom")
    )


    if (cluster_order == "size") {
        ## sorts clusters from big to small (# of cells per cluster)
        clus_sort_sizes <- sort(table(cell_metadata[[cluster_column]]))
        clus_sort_names <- names(clus_sort_sizes)
    } else if (cluster_order == "correlation") {
        ## sorts clusters based on their correlation
        subset_matrix <- create_cluster_matrix(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = feat_type,
            cluster_column = cluster_column,
            feat_subset = detected_feats,
            expression_values = values
        )

        cormatrix <- cor_flex(x = subset_matrix, method = cor_method)
        cordist <- stats::as.dist(1 - cormatrix, diag = TRUE, upper = TRUE)
        corclus <- stats::hclust(d = cordist, method = hclust_method)
        clus_names <- rownames(cormatrix)
        names(clus_names) <- seq_len(length(clus_names))
        clus_sort_names <- clus_names[corclus$order]
        clus_sort_names <- gsub(
            pattern = "cluster_",
            replacement = "", x = clus_sort_names
        )
    } else if (cluster_order == "custom") {
        ## sorts based on a given custom order
        if (is.null(cluster_custom_order)) {
            stop("\n custom order parameter is not given \n")
        }
        clus_sort_names <- cluster_custom_order
    }
    return(clus_sort_names)
}




# ggplot helper ####

#' @title aes_string2
#' @description
#' Internal replacement for the removed `ggplot2::aes_string()` pre-v4.0
#'
#' @param x,y,... List of name-value pairs in the form `aesthetic = variable`
#' describing which variables in the layer data should be mapped to which
#' aesthetics used by the paired geom/stat.
#'
#' * Mapped values should be provided as a single character input.
#' * Non-character inputs will not be altered and will be processed as usual.
#' * Expressions to be evaluated may be passed as a character string to be
#' interpreted via `str2lang()`
#'
#' @returns An S7 object representing a list with class mapping.
#' @examples
#' ggplot2::ggplot(iris) +
#' ggplot2::geom_point(aes_string2(x = "Sepal.Length", y = "Petal.Length"))
#' @export
aes_string2 <- function(x, y, ...) {
    mapping <- list(...)
    if (!missing(x))
        mapping["x"] <- list(x)
    if (!missing(y))
        mapping["y"] <- list(y)
    mapping <- lapply(mapping, function(x) {
        if (is.character(x)) {
            if (length(x) != 1L) {
                stop("[aes_string2] incompatible with len >1 character inputs",
                    call. = FALSE)
            }
            if (identical(x, "")) {
                stop("[aes_string2] empty string not allowed", call. = FALSE)
            }
            parsed <- str2lang(x)
            x <- if (is.symbol(parsed) || is.call(parsed)) {
                parsed
            }  else {
                as.name(x)
            }
        }
        x
    })
    do.call(ggplot2::aes, mapping)
}

# handle ggplot inputs for functions that may either
# append additional information to a ggplot object or be where the ggobject is
# first made.
# `ggobject` ggplot object or NULL
# returns ggplot object
.gg_input <- function(ggobject) {
    if (is.null(ggobject)) {
        return(ggplot2::ggplot())
    } else {
        checkmate::assert_class(ggobject, "ggplot")
        return(ggobject)
    }
}


# giotto point plotting ####

#' @title giotto_point
#' @name giotto_point
#' @description Plot a point scatter layer via one of ggplot::geom_point,
#' scattermore::geom_scattermore or scattermore::geom_scattermost
#' @param plot_method which scatter plotting method to use.
#' @param size point size parameter
#' @param ext SpatVector. Spatial extent of the plotting region. Used when
#' considering raster size for rasterization with scattermore and scattermost
#' @param scattermost_xy 2-column object with data. Used with scattermost since
#' it doesn't talk to ggplot
#' @param scattermost_color Color vector (or a single color). Used with
#' scattermost since it doesn't talk to ggplot
#' @param \dots geom_point parameters
#' @keywords internal
#' @noRd
#' @returns ggplot2::geom_point layer
giotto_point <- function(plot_method = c(
    "ggplot", "scattermore", "scattermost"),
    size = 1,
    ext,
    scattermost_xy = NULL,
    scattermost_color = NULL,
    ...) {
    plot_method <- match.arg(
        arg = plot_method,
        choices = c("ggplot", "scattermore", "scattermost")
    )

    px <- getOption("giotto.plot_point_raster", 5e5)
    if (!missing(ext)) {
        edim <- setNames(range(ext), NULL)
        ext_area <- prod(edim)
        px_per_unit <- px / ext_area
        pixels <- ceiling(edim * sqrt(px_per_unit))
    } else {
        pixels <- rep(ceiling(sqrt(px)), 2L)
    }

    # assemble argslist
    a <- list(...)

    if (plot_method %in% c("scattermore", "scattermost")) {
        a$pointsize <- size
        a$pixels <- pixels
        vmsg(.is_debug = TRUE, "rasterization dims:", pixels)
    } else {
        a$size <- size
    }

    # dispatch to plot method
    switch(plot_method,
        "ggplot" = {
            do.call(ggplot2::geom_point, a)
        },
        "scattermore" = {
            package_check(pkg_name = "scattermore", repository = "CRAN")
            do.call(scattermore::geom_scattermore, a)
        },
        "scattermost" = {
            package_check(pkg_name = "scattermore", repository = "CRAN")
            a$xy <- scattermost_xy
            a$color <- scattermost_color
            do.call(scattermore::geom_scattermost, a)
        }
    )
}


#' @name giotto_point_3d
#' @param pl plotly plot object
#' @param data dataframe-like object with data to plot
#' @param cell_color column in pl to plot. May be passed as NULL in which case
#' spots will default to `lightblue` in color.
#' @param color_as_factor logical. Whether values should be plotted as
#' categorical. Default is TRUE.
#' @param cell_color_code specific set of color hex codes to use in plotting
#' @param cell_color_gradient character. name of gradient to use
#' @param gradient_limits numeric of length 2. Numerical min and max values to
#' display in gradient scale.
#' @param point_size numeric. Size of points
#' @param data_other dataframe-like object with 'other' data points to plot
#' @param select_cells character. specific cell_IDs to select. Unselected will
#' be treated as 'other'
#' @param show_other_cells logical. whether to plot the 'other' cells.
#' @param other_cell_color character. Color code(s) to apply to 'other' cells.
#' @param other_point_size numeric. Size of points for 'other' cells
#' @param other_cell_alpha numeric. Alpha of 'other' cells
#' @param instrs giottoInstructions
#' @keywords internal
#' @noRd
giotto_point_3d <- function(pl,
    data,
    cell_color = NULL,
    color_as_factor = TRUE,
    cell_color_code = NULL,
    cell_color_gradient = NULL,
    gradient_limits = NULL,
    gradient_style = "divergent",
    gradient_midpoint = NULL,
    point_size = 3,
    point_alpha = 1,
    data_other = NULL,
    select_cells = NULL,
    show_other_cells = TRUE,
    other_cell_color = "lightgrey",
    other_point_size = 0.5,
    other_cell_alpha = 3,
    instrs) {
    # plotly params list init & static params
    # ** toplevel ** #
    trace_params <- trace_params_other <- list(
        type = "scatter3d",
        mode = "markers",
        x = ~sdimx,
        y = ~sdimy,
        z = ~sdimz
    )
    trace_params$colors <- "lightblue" # default with `cell_color` = NULL
    trace_params$opacity <- point_alpha
    trace_params_other$name <- "unselected cells"
    trace_params_other$opacity <- other_cell_alpha

    # ** point aes ** #
    marker_params <- marker_params_other <- list()
    marker_params$size <- point_size
    marker_params_other$size <- other_point_size
    marker_params_other$color <- other_cell_color

    # finalize data and color
    if (!is.null(cell_color)) {
        if (!cell_color %in% colnames(data)) {
            message(sprintf("`cell_color` '%s' does not exist!", cell_color))
            return(pl) # return early
        }

        # finalize color scale
        if (is.null(cell_color_code)) {
            if (color_as_factor) { # categorical
                number_colors <- length(unique(
                    data[[cell_color]]
                ))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instrs
                )(n = number_colors)
            } else { # continuous
                cell_color_code <- set_default_color_continuous_cell(
                    instrs = instrs,
                    style = gradient_style,
                    midpoint = gradient_midpoint,
                    colors = cell_color_gradient
                )$palette(seq(0, 1, length.out = 100L))
            }
        }

        # finalize data
        if (color_as_factor) { # categorical
            data[[cell_color]] <- as.factor(data[[cell_color]])
        } else { # continuous
            # assign gradient limits if needed
            if (!is.null(gradient_limits)) {
                checkmate::assert_numeric(gradient_limits, len = 2L)
                lower_lim <- gradient_limits[[1L]]
                upper_lim <- gradient_limits[[2L]]
                data[, (cell_color) :=
                    scales::oob_squish(get(cell_color), gradient_limits)]
            }
        }
        # apply non-default color settings
        trace_params$color <- data[[cell_color]]
        trace_params$colors <- cell_color_code
    }
    # apply finalized data and marker
    trace_params$data <- data
    trace_params$marker <- marker_params

    # other points
    if (!is.null(select_cells) && isTRUE(show_other_cells)) {
        trace_params$name <- "selected_cells"

        trace_params_other$data <- data_other
        trace_params_other$marker <- marker_params_other
        pl <- do.call(
            plotly::add_trace,
            args = c(list(p = pl), trace_params_other)
        )
    }

    pl <- do.call(plotly::add_trace, args = c(list(p = pl), trace_params))

    return(pl)
}




# rescale values ####

# based on ggplot2 internal
mid_rescaler <- function(mid) {
    function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
        scales::rescale_mid(x, to, from, mid)
    }
}






# plotly helper ####

#' @title plotly_network
#' @name plotly_network
#' @description provide network segment to draw in 3D plot_ly()
#' @param network network object
#' @param x default to "sdimx_begin"
#' @param y default to "sdimy_begin"
#' @param z default to "sdimz_begin"
#' @param x_end default to "sdimx_end"
#' @param y_end default to "sdimy_end"
#' @param z_end default to "sdimz_end"
#'
#' @returns edges in network as data.table
#'
#' @export
plotly_network <- function(
        network,
        x = "sdimx_begin",
        y = "sdimy_begin",
        z = "sdimz_begin",
        x_end = "sdimx_end",
        y_end = "sdimy_end",
        z_end = "sdimz_end") {
    edges <- data.table::data.table(
        edge_id = seq_len(3 * dim(network)[1]),
        x = 0,
        y = 0,
        z = 0
    )

    edges[edges$edge_id %% 3 == 1]$x <- as.double(network[[x]])
    edges[edges$edge_id %% 3 == 1]$y <- as.double(network[[y]])
    edges[edges$edge_id %% 3 == 1]$z <- as.double(network[[z]])

    edges[edges$edge_id %% 3 == 2]$x <- as.double(network[[x_end]])
    edges[edges$edge_id %% 3 == 2]$y <- as.double(network[[y_end]])
    edges[edges$edge_id %% 3 == 2]$z <- as.double(network[[z_end]])

    edges[edges$edge_id %% 3 == 0]$x <- NA
    edges[edges$edge_id %% 3 == 0]$y <- NA
    edges[edges$edge_id %% 3 == 0]$z <- NA

    return(edges)
}


#' @title plotly_grid
#' @name plotly_grid
#' @description provide grid segment to draw in plot_ly()
#'
#' @param x_start default to "x_start"
#' @param y_start default to "y_start"
#' @param x_end default to "x_end"
#' @param y_end default to "y_end"
#' @param spatial_grid spatial_grid in giotto object
#'
#' @returns edges in spatial grid as data.table()
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#'
#' library(GiottoClass)
#' g <- createSpatialGrid(
#'     gobject = g,
#'     sdimx_stepsize = 400,
#'     sdimy_stepsize = 400,
#'     minimum_padding = 0
#' )
#'
#' my_spatial_grid <- getSpatialGrid(g)
#'
#' plotly_grid(my_spatial_grid)
#'
#' @export
plotly_grid <- function(spatial_grid,
    x_start = "x_start",
    y_start = "y_start",
    x_end = "x_end",
    y_end = "y_end") {
    edge_num <- length(unique(spatial_grid[[x_start]])) +
        length(unique(spatial_grid[[y_start]])) + 2
    x_line <- unique(as.numeric(unlist(spatial_grid[, c(x_start, x_end)])))
    y_line <- unique(as.numeric(unlist(spatial_grid[, c(y_start, y_end)])))

    x_min <- min(spatial_grid[[x_start]])
    x_max <- max(spatial_grid[[x_end]])

    y_min <- min(spatial_grid[[y_start]])
    y_max <- max(spatial_grid[[y_end]])

    edges <- data.table::data.table(
        edge_id = seq_len(edge_num),
        x = 0, y = 0, x_end = 0, y_end = 0
    )

    edges[seq_len(length(x_line)), ]$x <- x_line
    edges[seq_len(length(x_line)), ]$x_end <- x_line
    edges[seq_len(length(x_line)), ]$y <- y_min
    edges[seq_len(length(x_line)), ]$y_end <- y_max

    edges[(length(x_line) + 1):edge_num, ]$x <- x_min
    edges[(length(x_line) + 1):edge_num, ]$x_end <- x_max
    edges[(length(x_line) + 1):edge_num, ]$y <- y_line
    edges[(length(x_line) + 1):edge_num, ]$y_end <- y_line

    return(edges)
}


#' @title plotly_axis_scale_3D
#' @name plotly_axis_scale_3D
#' @description adjust the axis scale in 3D plotly plot
#' @param cell_locations spatial_loc in giotto object
#' @param sdimx x axis of cell spatial location
#' @param sdimy y axis of cell spatial location
#' @param sdimz z axis of cell spatial location
#' @param mode axis adjustment mode
#' @param custom_ratio set the ratio artificially
#' @returns edges in spatial grid as data.table()
#' @examples
#' my_cell_locations <- data.frame(
#'     x = sample(10), y = sample(10),
#'     z = sample(10)
#' )
#' plotly_axis_scale_3D(my_cell_locations)
#'
#' @export
plotly_axis_scale_3D <- function(cell_locations,
    sdimx = NULL,
    sdimy = NULL,
    sdimz = NULL,
    mode = c("cube", "real", "custom"),
    custom_ratio = NULL) {
    mode <- match.arg(mode, c("cube", "real", "custom"))
    if (mode == "real") {
        x_ratio <- max(cell_locations[[sdimx]]) - min(cell_locations[[sdimx]])
        y_ratio <- max(cell_locations[[sdimy]]) - min(cell_locations[[sdimy]])
        z_ratio <- max(cell_locations[[sdimz]]) - min(cell_locations[[sdimz]])

        min_size <- min(x_ratio, y_ratio, z_ratio)
        x_ratio <- round(x_ratio / min_size)
        y_ratio <- round(y_ratio / min_size)
        z_ratio <- round(z_ratio / min_size)
    } else if (mode == "cube") {
        x_ratio <- 1
        y_ratio <- 1
        z_ratio <- 1
    } else {
        if (is.null(custom_ratio) | length(custom_ratio) < 3) {
            stop("\n Please specify your costom axis scale,
                or choose axis_scale = \"real\"/\"cube\"\n")
        } else {
            x_ratio <- custom_ratio[1]
            y_ratio <- custom_ratio[2]
            z_ratio <- custom_ratio[3]
        }
    }
    ratio <- list(x_ratio, y_ratio, z_ratio)
    return(ratio)
}


#' @title plotly_axis_scale_2D
#' @name plotly_axis_scale_2D
#' @description adjust the axis scale in 2D plotly plot
#' @param cell_locations spatial_loc in giotto object
#' @param sdimx x axis of cell spatial location
#' @param sdimy y axis of cell spatial location
#' @param mode axis adjustment mode
#' @param custom_ratio set the ratio artificially
#' @returns edges in spatial grid as data.table()
#' @examples
#' my_cell_locations <- data.frame(x = sample(10), y = sample(10))
#' plotly_axis_scale_2D(my_cell_locations)
#'
#' @export
plotly_axis_scale_2D <- function(cell_locations,
    sdimx = NULL,
    sdimy = NULL,
    mode = c("cube", "real", "custom"),
    custom_ratio = NULL) {
    mode <- match.arg(mode, c("cube", "real", "custom"))
    if (mode == "real") {
        x_ratio <- max(cell_locations[[sdimx]]) - min(cell_locations[[sdimx]])
        y_ratio <- max(cell_locations[[sdimy]]) - min(cell_locations[[sdimy]])

        min_size <- min(x_ratio, y_ratio)
        x_ratio <- round(x_ratio / min_size)
        y_ratio <- round(y_ratio / min_size)
    } else if (mode == "cube") {
        x_ratio <- 1
        y_ratio <- 1
    } else {
        if (is.null(custom_ratio) | length(custom_ratio) < 2) {
            stop("\n Please specify your costom axis scale,
                or choose axis_scale = \"real\"/\"cube\"\n")
        } else {
            x_ratio <- custom_ratio[1]
            y_ratio <- custom_ratio[2]
        }
    }
    ratio <- list(x_ratio, y_ratio)
    return(ratio)
}
