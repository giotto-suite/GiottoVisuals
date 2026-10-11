# --------------------------------------------------------------------------- #
# INTERNAL                                                                    #
# ggplot2 helper functions that are used to create layers of information that #
# are used in the final plot.                                                 #
# --------------------------------------------------------------------------- #

# TODO clean up this code a little more
#   (there should only be one geom_point() call)
# TODO support dimFeatPlot()

# spatial ####

#' @title plot_spat_point_layer_ggplot
#' @name plot_spat_point_layer_ggplot
#' @description create ggplot point layer for spatial coordinates
#' @param ggobject ggplot object
#' @inheritParams plot_params
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_locations_metadata_selected annotated location from
#' selected cells
#' @param cell_locations_metadata_other annotated location from
#' non-selected cells
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @noRd
plot_spat_point_layer_ggplot <- function(
        ggobject,
        instrs = NULL,
        sdimx = NULL,
        sdimy = NULL,
        plot_method = c("ggplot", "scattermore"), # does not work for fill types
        cell_locations_metadata_selected,
        cell_locations_metadata_other,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = "divergent",
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        point_size = 2,
        point_alpha = 1,
        point_border_col = "lightgrey",
        point_border_stroke = 0.1,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        show_legend = TRUE,
        ...) {
    ## specify spatial dimensions first
    if (is.null(sdimx) || is.null(sdimy)) {
        warning(wrap_txt("plot_method = ggplot,
                        but spatial dimensions for sdimx and/or sdimy
                        are not specified. It will default to the
                        'sdimx' and 'sdimy' "))
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }

    plot_method <- match.arg(plot_method, c("ggplot", "scattermore"))

    ### point parameters ##
    point_size <- as.numeric(point_size)
    point_alpha <- as.numeric(point_alpha)
    point_border_stroke <- as.numeric(point_border_stroke)

    center_point_size <- as.numeric(center_point_size)
    center_point_border_stroke <- as.numeric(center_point_border_stroke)

    label_size <- as.numeric(label_size)
    other_point_size <- as.numeric(other_point_size)


    ## ggplot object
    pl <- ggobject

    ## first plot other non-selected cells
    if ((!is.null(select_cells) || !is.null(select_cell_groups)) &&
        isTRUE(show_other_cells)) {
        pl <- pl + ggplot2::geom_point(
            data = cell_locations_metadata_other,
            mapping = aes_string2(x = sdimx, sdimy),
            color = other_cell_color,
            show.legend = FALSE,
            size = other_point_size,
            alpha = point_alpha
        )
    }


    ## order of color
    # 1. if NULL then default to lightblue
    # 2. if character vector
    # 2.1 if length of cell_color is longer than 1 and has colors
    # 2.2 if not part of metadata then suppose its color
    # 2.3 part of metadata
    # 2.3.1 numerical column
    # 2.3.2 factor column or character to factor

    # cell color default
    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl + giotto_point(
            plot_method = plot_method,
            data = cell_locations_metadata_selected,
            mapping = aes_string2(x = sdimx, y = sdimy),
            show.legend = show_legend,
            shape = 21,
            fill = cell_color,
            size = point_size,
            stroke = point_border_stroke,
            color = point_border_col,
            alpha = point_alpha
        )
    } else if (length(cell_color) > 1L) {
        if (is.numeric(cell_color) || is.factor(cell_color)) {
            # cell_color is numeric
            if (nrow(cell_locations_metadata_selected) != length(cell_color)) {
                stop("\n vector needs to be the same lengths as
                    number of cells \n")
            }
            cell_locations_metadata_selected[["temp_color"]] <- cell_color

            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = cell_locations_metadata_selected,
                mapping = aes_string2(
                    x = sdimx,
                    y = sdimy,
                    fill = "temp_color"),
                show.legend = show_legend,
                shape = 21,
                size = point_size,
                color = point_border_col,
                stroke = point_border_stroke,
                alpha = point_alpha
            )
        } else if (is.character(cell_color)) { # cell_color is hex codes
            if (!all(.is_color_code(cell_color))) {
                stop("cell_color is not numeric,
                    a factor or vector of colors \n")
            }

            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = cell_locations_metadata_selected,
                mapping = aes_string2(x = sdimx, y = sdimy),
                show.legend = show_legend,
                shape = 21,
                fill = cell_color,
                size = point_size,
                color = point_border_col,
                stroke = point_border_stroke,
                alpha = point_alpha
            )
        }
    } else if (is.character(cell_color)) {
        if (!cell_color %in% colnames(cell_locations_metadata_selected)) {
            if (!.is_color_code(cell_color)) {
                stop(cell_color, " is not a color or a column name \n")
            }
            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = cell_locations_metadata_selected,
                mapping = aes_string2(x = sdimx, y = sdimy),
                show.legend = show_legend,
                shape = 21,
                fill = cell_color,
                size = point_size,
                color = point_border_col,
                stroke = point_border_stroke,
                alpha = point_alpha
            )
        } else {
            class_cell_color <-
                class(cell_locations_metadata_selected[[cell_color]])

            if ((class_cell_color %in% c("integer", "numeric")) &&
                !isTRUE(color_as_factor)) {
                # set upper and lower limits
                if (!is.null(gradient_limits) &&
                    is.vector(gradient_limits) &&
                    length(gradient_limits) == 2) {
                    lower_lim <- gradient_limits[[1]]
                    upper_lim <- gradient_limits[[2]]

                    numeric_data <-
                        cell_locations_metadata_selected[[cell_color]]
                    limit_numeric_data <-
                        ifelse(numeric_data > upper_lim, upper_lim,
                            ifelse(numeric_data < lower_lim, lower_lim,
                                numeric_data
                            )
                        )
                    cell_locations_metadata_selected[[cell_color]] <-
                        limit_numeric_data
                }

                pl <- pl + giotto_point(
                    plot_method = plot_method,
                    data = cell_locations_metadata_selected,
                    mapping = aes_string2(
                        x = sdimx,
                        y = sdimy,
                        fill = cell_color),
                    show.legend = show_legend,
                    shape = 21,
                    size = point_size,
                    color = point_border_col,
                    stroke = point_border_stroke,
                    alpha = point_alpha
                )
            } else {
                # convert character or numeric to factor
                if (isTRUE(color_as_factor)) {
                    factor_data <-
                        factor(cell_locations_metadata_selected[[cell_color]])
                    cell_locations_metadata_selected[[cell_color]] <-
                        factor_data
                }

                # if you want to show centers or labels then calculate centers
                if (isTRUE(show_cluster_center) || isTRUE(show_center_label)) {
                    annotated_DT_centers <-
                        cell_locations_metadata_selected[, .(
                            center_1 = stats::median(get("sdimx")),
                            center_2 = stats::median(get("sdimy"))
                        ), by = cell_color]
                    factor_center_data <-
                        factor(annotated_DT_centers[[cell_color]])
                    annotated_DT_centers[[cell_color]] <- factor_center_data
                }

                pl <- pl + giotto_point(
                    plot_method = plot_method,
                    data = cell_locations_metadata_selected,
                    mapping = aes_string2(
                        x = sdimx,
                        y = sdimy,
                        fill = cell_color),
                    show.legend = show_legend,
                    shape = 21,
                    size = point_size,
                    color = point_border_col,
                    stroke = point_border_stroke,
                    alpha = point_alpha
                )


                ## plot centers
                if (isTRUE(show_cluster_center) &&
                    (isTRUE(color_as_factor) ||
                        class_cell_color %in% c("character", "factor"))) {
                    pl <- pl + giotto_point(
                        plot_method = plot_method,
                        data = annotated_DT_centers,
                        mapping = aes_string2(
                            x = "center_1",
                            y = "center_2",
                            fill = cell_color
                        ),
                        color = center_point_border_col,
                        stroke = center_point_border_stroke,
                        size = center_point_size,
                        shape = 21,
                        alpha = point_alpha
                    )
                }

                ## plot labels
                if (isTRUE(show_center_label)) {
                    pl <- pl + ggrepel::geom_text_repel(
                        data = annotated_DT_centers,
                        mapping = aes_string2(
                            x = "center_1",
                            y = "center_2",
                            label = cell_color
                        ),
                        size = label_size,
                        fontface = label_fontface
                    )
                }
            }

            ## specify colors to use
            if (!is.null(cell_color_code)) {
                pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (isTRUE(color_as_factor)) {
                number_colors <- length(unique(factor_data))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instrs
                )(n = number_colors)
                names(cell_color_code) <- unique(factor_data)
                pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (!isTRUE(color_as_factor)) {
                if (is.null(gradient_midpoint)) {
                    gradient_midpoint <- stats::median(
                        cell_locations_metadata_selected[[cell_color]]
                    )
                }

                pl <- pl + set_default_color_continuous_cell(
                    colors = cell_color_gradient,
                    instrs = instrs,
                    style = gradient_style,
                    midpoint = gradient_midpoint,
                    type = "fill",
                )
            }
        }
    }
    return(pl)
}


#' @title plot_spat_point_layer_ggplot_noFILL
#' @name plot_spat_point_layer_ggplot_noFILL
#' @description create ggplot point layer for spatial coordinates without
#' borders
#' @param ggobject ggplot object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_locations_metadata_selected annotated location from selected
#' cells
#' @param cell_locations_metadata_other annotated location from non-selected
#' cells
#' @inheritParams plot_params
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @noRd
plot_spat_point_layer_ggplot_noFILL <- function(
        ggobject,
        instrs = NULL,
        sdimx = NULL,
        sdimy = NULL,
        plot_method = c("ggplot", "scattermore"),
        cell_locations_metadata_selected,
        cell_locations_metadata_other,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = "divergent",
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        point_size = 2,
        point_alpha = 1,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        label_size = 4,
        label_fontface = "bold",
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        show_legend = TRUE,
        ...) {
    ## specify spatial dimensions first
    if (is.null(sdimx) || is.null(sdimy)) {
        warning(wrap_txt("plot_method = ggplot, but spatial dimensions
                        for sdimx and/or sdimy are not specified.
                        It will default to the 'sdimx' and 'sdimy'"))
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }

    plot_method <- match.arg(plot_method, c("ggplot", "scattermore"))

    ## point parameters ##
    point_size <- as.numeric(point_size)
    point_alpha <- as.numeric(point_alpha)
    center_point_size <- as.numeric(center_point_size)
    label_size <- as.numeric(label_size)
    other_point_size <- as.numeric(other_point_size)


    ## ggplot object
    pl <- ggobject

    ## first plot other non-selected cells
    if ((!is.null(select_cells) || !is.null(select_cell_groups)) &&
        isTRUE(show_other_cells)) {
        pl <- pl + ggplot2::geom_point(
            data = cell_locations_metadata_other,
            aes_string2(x = sdimx, sdimy),
            color = other_cell_color,
            show.legend = FALSE,
            size = other_point_size,
            alpha = point_alpha
        )
    }


    ## order of color
    # 1. if NULL then default to lightblue
    # 2. if character vector
    # 2.1 if length of cell_color is longer than 1 and has colors
    # 2.2 if not part of metadata then suppose its color
    # 2.3 part of metadata
    # 2.3.1 numerical column
    # 2.3.2 factor column or character to factor

    # cell color default
    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl + giotto_point(
            plot_method = plot_method,
            data = cell_locations_metadata_selected,
            mapping = aes_string2(x = sdimx, y = sdimy),
            show.legend = show_legend,
            shape = 19,
            color = cell_color,
            size = point_size,
            alpha = point_alpha,
            ...
        )
    } else if (length(cell_color) > 1) {
        if (is.numeric(cell_color) || is.factor(cell_color)) {
            if (nrow(cell_locations_metadata_selected) != length(cell_color)) {
                stop("\n vector needs to be the same lengths as number
                    of cells \n")
            }
            cell_locations_metadata_selected[["temp_color"]] <- cell_color

            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = cell_locations_metadata_selected,
                mapping = aes_string2(
                    x = sdimx,
                    y = sdimy,
                    color = "temp_color"),
                show.legend = show_legend,
                shape = 19,
                size = point_size,
                alpha = point_alpha,
                ...
            )
        } else if (is.character(cell_color)) {
            if (!all(.is_color_code(cell_color))) {
                stop("cell_color is not numeric, a factor or vector of
                    colors \n")
            }

            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = cell_locations_metadata_selected,
                mapping = aes_string2(x = sdimx, y = sdimy),
                show.legend = show_legend,
                shape = 19,
                color = cell_color,
                size = point_size,
                alpha = point_alpha,
                ...
            )
        }
    } else if (is.character(cell_color)) {
        if (!cell_color %in% colnames(cell_locations_metadata_selected)) {
            if (!.is_color_code(cell_color)) {
                stop(cell_color, " is not a color or a column name \n")
            }
            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = cell_locations_metadata_selected,
                mapping = aes_string2(x = sdimx, y = sdimy),
                show.legend = show_legend,
                shape = 19,
                color = cell_color,
                size = point_size,
                alpha = point_alpha,
                ...
            )
        } else {
            class_cell_color <- class(
                cell_locations_metadata_selected[[cell_color]]
            )

            if (class_cell_color %in% c("integer", "numeric") &&
                !isTRUE(color_as_factor)) {
                # set upper and lower limits
                if (!is.null(gradient_limits) &&
                    is.vector(gradient_limits) &&
                    length(gradient_limits) == 2) {
                    lower_lim <- gradient_limits[[1]]
                    upper_lim <- gradient_limits[[2]]

                    numeric_data <-
                        cell_locations_metadata_selected[[cell_color]]
                    limit_numeric_data <- ifelse(
                        numeric_data > upper_lim, upper_lim,
                        ifelse(numeric_data < lower_lim, lower_lim,
                            numeric_data
                        )
                    )
                    cell_locations_metadata_selected[[cell_color]] <-
                        limit_numeric_data
                }

                pl <- pl + giotto_point(
                    plot_method = plot_method,
                    data = cell_locations_metadata_selected,
                    mapping = aes_string2(
                        x = sdimx,
                        y = sdimy,
                        color = cell_color),
                    show.legend = show_legend,
                    shape = 19,
                    size = point_size,
                    alpha = point_alpha,
                    ...
                )
            } else {
                # convert character or numeric to factor
                if (isTRUE(color_as_factor)) {
                    factor_data <- factor(
                        cell_locations_metadata_selected[[cell_color]]
                    )
                    cell_locations_metadata_selected[[cell_color]] <-
                        factor_data
                }

                # if you want to show centers or labels then calculate centers
                if (isTRUE(show_cluster_center) || isTRUE(show_center_label)) {
                    annotated_DT_centers <-
                        cell_locations_metadata_selected[, .(
                            center_1 = stats::median(get("sdimx")),
                            center_2 = stats::median(get("sdimy"))
                        ), by = cell_color]
                    factor_center_data <-
                        factor(annotated_DT_centers[[cell_color]])
                    annotated_DT_centers[[cell_color]] <- factor_center_data
                }

                pl <- pl + giotto_point(
                    plot_method = plot_method,
                    data = cell_locations_metadata_selected,
                    mapping = aes_string2(
                        x = sdimx,
                        y = sdimy,
                        color = cell_color),
                    show.legend = show_legend,
                    shape = 19,
                    size = point_size,
                    alpha = point_alpha,
                    ...
                )


                ## plot centers
                if (isTRUE(show_cluster_center) &&
                    (isTRUE(color_as_factor) ||
                        class_cell_color %in% c("character", "factor"))) {
                    pl <- pl + giotto_point(
                        plot_method = plot_method,
                        data = annotated_DT_centers,
                        mapping = aes_string2(
                            x = "center_1",
                            y = "center_2",
                            color = cell_color
                        ),
                        size = center_point_size,
                        shape = 19,
                        alpha = point_alpha,
                        ...
                    )
                }

                ## plot labels
                if (isTRUE(show_center_label)) {
                    pl <- pl + ggrepel::geom_text_repel(
                        data = annotated_DT_centers,
                        aes_string2(
                            x = "center_1",
                            y = "center_2",
                            label = cell_color
                        ),
                        size = label_size,
                        fontface = label_fontface,
                        alpha = point_alpha
                    )
                }
            }

            ## specificy colors to use
            if (!is.null(cell_color_code)) {
                pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)
            } else if (isTRUE(color_as_factor)) {
                number_colors <- length(unique(factor_data))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instrs
                )(n = number_colors)
                names(cell_color_code) <- unique(factor_data)
                pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)
            } else if (!isTRUE(color_as_factor)) {
                if (is.null(gradient_midpoint)) {
                    gradient_midpoint <- stats::median(
                        cell_locations_metadata_selected[[cell_color]]
                    )
                }

                pl <- pl + set_default_color_continuous_cell(
                    colors = cell_color_gradient,
                    instrs = instrs,
                    style = gradient_style,
                    midpoint = gradient_midpoint,
                    type = "color",
                )
            }
        }
    }
    return(pl)
}



#' @title plot_spat_voronoi_layer_ggplot
#' @name plot_spat_voronoi_layer_ggplot
#' @description creat ggplot point layer for spatial coordinates without borders
#' @param ggobject ggplot object
#' @inheritParams plot_params
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_locations_metadata_selected annotated location from
#' selected cells
#' @param cell_locations_metadata_other annotated location from
#' non-selected cells
#' @param vor_border_color borde colorr of voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @noRd
plot_spat_voronoi_layer_ggplot <- function(
        ggobject,
        instrs = NULL,
        sdimx = NULL,
        sdimy = NULL,
        cell_locations_metadata_selected,
        cell_locations_metadata_other,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = NULL,
        gradient_style = "divergent",
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        point_size = 2,
        point_alpha = 1,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        label_size = 4,
        label_fontface = "bold",
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 1,
        background_color = "white",
        vor_border_color = "white",
        vor_max_radius = 200,
        vor_alpha = 1,
        show_legend = TRUE,
        ...) {
    ## specify spatial dimensions first
    if (is.null(sdimx) | is.null(sdimy)) {
        warning("plot_method = ggplot, but spatial dimensions for sdimx
                and/or sdimy are not specified. \n
                It will default to the 'sdimx' and 'sdimy' ")
        sdimx <- "sdimx"
        sdimy <- "sdimy"
    }

    ## ggplot object
    pl <- ggobject



    ## order of color
    # 1. if NULL then default to lightblue
    # 2. if character vector
    # 2.1 if length of cell_color is longer than 1 and has colors
    # 2.2 if not part of metadata then suppose its color
    # 2.3 part of metadata
    # 2.3.1 numerical column
    # 2.3.2 factor column or character to factor

    # data.table variables
    temp_color <- NULL

    # cell color default
    if (is.null(cell_color)) {
        ## 1. default colors when no colors are assigned ##

        cell_color <- "lightblue"
        cell_locations_metadata_selected[, "temp_color" := "selected"]

        if (!is.null(cell_locations_metadata_other)) {
            cell_locations_metadata_other[, "temp_color" := "other"]
        }


        combn_cell_locations_metadata <-
            rbind(
                cell_locations_metadata_selected,
                cell_locations_metadata_other
            )

        pl <- pl + ggforce::geom_voronoi_tile(
            data = combn_cell_locations_metadata,
            aes(
                x = sdimx,
                y = sdimy,
                group = -1L,
                fill = as.factor(temp_color)
            ),
            colour = vor_border_color,
            max.radius = vor_max_radius,
            show.legend = show_legend,
            alpha = vor_alpha
        )

        if (show_other_cells == TRUE) {
            pl <- pl +
                ggplot2::scale_fill_manual(values = c(
                    selected = cell_color,
                    other = other_cell_color
                ))
        } else {
            pl <- pl +
                ggplot2::scale_fill_manual(values = c(
                    selected = cell_color,
                    other = background_color
                ))
        }

        # theme specific changes
        pl <- pl + theme(legend.title = element_blank())
    } else if (length(cell_color) > 1L) {
        ## 2. continuous vector to convert to colors ##
        if (is.numeric(cell_color) || is.factor(cell_color)) {
            if (nrow(cell_locations_metadata_selected) != length(cell_color)) {
                stop("\n vector needs to be the same lengths as number
                    of cells \n")
            }

            cell_locations_metadata_selected[["temp_color"]] <- cell_color
            if (!is.null(cell_locations_metadata_other)) {
                cell_locations_metadata_other[["temp_color"]] <- NA
            }
            combn_cell_locations_metadata <-
                rbind(
                    cell_locations_metadata_selected,
                    cell_locations_metadata_other
                )

            pl <- pl + ggforce::geom_voronoi_tile(
                data = combn_cell_locations_metadata,
                aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                colour = vor_border_color, max.radius = vor_max_radius,
                show.legend = show_legend,
                alpha = vor_alpha
            )

            if (is.null(gradient_midpoint)) {
                gradient_midpoint <- stats::median(
                    cell_locations_metadata_selected[["temp_color"]]
                )
            }

            mybg_color <- ifelse(show_other_cells == TRUE, other_cell_color,
                background_color
            )

            pl <- pl + set_default_color_continuous_cell(
                colors = cell_color_gradient,
                instrs = instrs,
                midpoint = gradient_midpoint,
                style = gradient_style,
                na.value = mybg_color,
                type = "fill"
            )

            # theme specific changes
            pl <- pl + theme(legend.title = element_blank())
        } else if (is.character(cell_color)) {
            ## 3. character vector to convert to colors ##

            if (!all(.is_color_code(cell_color))) {
                stop("cell_color is not numeric, a factor or vector of
                    colors \n")
            }

            if (nrow(cell_locations_metadata_selected) != length(cell_color)) {
                stop("\n vector needs to be the same lengths as number
                    of cells \n")
            }

            other_cell_color <- ifelse(isTRUE(show_other_cells),
                other_cell_color, background_color
            )

            cell_locations_metadata_selected[["temp_color"]] <- cell_color
            if (!is.null(cell_locations_metadata_other)) {
                cell_locations_metadata_other[["temp_color"]] <-
                    other_cell_color
            }
            combn_cell_locations_metadata <-
                rbind(
                    cell_locations_metadata_selected,
                    cell_locations_metadata_other
                )

            pl <- pl + ggforce::geom_voronoi_tile(
                data = combn_cell_locations_metadata,
                aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                colour = vor_border_color,
                max.radius = vor_max_radius,
                show.legend = show_legend,
                alpha = vor_alpha
            )

            my_color_code <-
                unique(combn_cell_locations_metadata[["temp_color"]])
            names(my_color_code) <- my_color_code

            pl <- pl + ggplot2::scale_fill_manual(values = my_color_code)

            # theme specific changes
            pl <- pl + theme(legend.title = element_blank())
        }
    } else if (is.character(cell_color)) {
        if (!cell_color %in% colnames(cell_locations_metadata_selected)) {
            if (!.is_color_code(cell_color)) {
                stop(cell_color, " is not a color or a column name \n")
            }

            ## 4. use a specific color ##
            other_cell_color <- ifelse(show_other_cells == TRUE,
                other_cell_color, background_color
            )

            cell_locations_metadata_selected[["temp_color"]] <- "selected"
            if (!is.null(cell_locations_metadata_other)) {
                cell_locations_metadata_other[["temp_color"]] <- "other"
            }
            combn_cell_locations_metadata <-
                rbind(
                    cell_locations_metadata_selected,
                    cell_locations_metadata_other
                )

            pl <- pl + ggforce::geom_voronoi_tile(
                data = combn_cell_locations_metadata,
                aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                colour = vor_border_color,
                max.radius = vor_max_radius,
                show.legend = show_legend,
                alpha = vor_alpha
            )

            my_color_code <- unique(
                combn_cell_locations_metadata[["temp_color"]]
            )
            names(my_color_code) <- my_color_code
            pl <- pl +
                ggplot2::scale_fill_manual(values = c(
                    selected = cell_color,
                    other = other_cell_color
                ))

            # theme specific changes
            pl <- pl + theme(legend.title = element_blank())
        } else {
            class_cell_color <-
                class(cell_locations_metadata_selected[[cell_color]])

            if ((class_cell_color %in% c("integer", "numeric")) &&
                !isTRUE(color_as_factor)) {
                ## 5. use continuous column from metadata ##

                # set upper and lower limits
                if (!is.null(gradient_limits) &&
                    is.vector(gradient_limits) &&
                    length(gradient_limits) == 2L) {
                    lower_lim <- gradient_limits[[1]]
                    upper_lim <- gradient_limits[[2]]

                    numeric_data <-
                        cell_locations_metadata_selected[[cell_color]]

                    limit_numeric_data <- ifelse(numeric_data > upper_lim,
                        upper_lim,
                        ifelse(numeric_data < lower_lim,
                            lower_lim, numeric_data
                        )
                    )
                    cell_locations_metadata_selected[[cell_color]] <-
                        limit_numeric_data
                }

                cell_locations_metadata_selected[["temp_color"]] <-
                    cell_locations_metadata_selected[[cell_color]]
                if (!is.null(cell_locations_metadata_other)) {
                    cell_locations_metadata_other[["temp_color"]] <- NA
                }
                combn_cell_locations_metadata <-
                    rbind(
                        cell_locations_metadata_selected,
                        cell_locations_metadata_other
                    )

                pl <- pl + ggforce::geom_voronoi_tile(
                    data = combn_cell_locations_metadata,
                    aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                    colour = vor_border_color,
                    max.radius = vor_max_radius,
                    show.legend = show_legend,
                    alpha = vor_alpha
                )

                mybg_color <- ifelse(isTRUE(show_other_cells),
                    other_cell_color, background_color
                )

                if (is.null(gradient_midpoint)) {
                    gradient_midpoint <- stats::median(
                        cell_locations_metadata_selected[["temp_color"]]
                    )
                }

                pl <- pl + set_default_color_continuous_cell(
                    colors = cell_color_gradient,
                    instrs = instrs,
                    midpoint = gradient_midpoint,
                    style = gradient_style,
                    na.value = mybg_color,
                    name = cell_color,
                    type = "fill"
                )
            } else {
                ## 6. use factor or character column from metadata ##
                # convert character or numeric to factor
                if (color_as_factor == TRUE) {
                    factor_data <- factor(
                        cell_locations_metadata_selected[[cell_color]]
                    )
                    cell_locations_metadata_selected[[cell_color]] <-
                        factor_data
                }

                # if you want to show centers or labels then calculate centers
                if (show_cluster_center == TRUE || show_center_label == TRUE) {
                    annotated_DT_centers <-
                        cell_locations_metadata_selected[, .(
                            center_1 = stats::median(get("sdimx")),
                            center_2 = stats::median(get("sdimy"))
                        ), by = cell_color]
                    factor_center_data <- factor(
                        annotated_DT_centers[[cell_color]]
                    )
                    annotated_DT_centers[[cell_color]] <- factor_center_data
                }

                cell_locations_metadata_selected[["temp_color"]] <-
                    cell_locations_metadata_selected[[cell_color]]
                if (!is.null(cell_locations_metadata_other)) {
                    cell_locations_metadata_other[["temp_color"]] <- "other"
                }
                combn_cell_locations_metadata <- rbind(
                    cell_locations_metadata_selected,
                    cell_locations_metadata_other
                )

                pl <- pl + ggforce::geom_voronoi_tile(
                    data = combn_cell_locations_metadata,
                    aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                    colour = vor_border_color,
                    max.radius = vor_max_radius,
                    show.legend = show_legend,
                    alpha = vor_alpha
                )


                other_cell_color <- ifelse(isTRUE(show_other_cells),
                    other_cell_color, background_color
                )

                ## specify colors to use
                if (!is.null(cell_color_code)) {
                    cell_color_code[["other"]] <- other_cell_color
                    pl <- pl + ggplot2::scale_fill_manual(
                        values = cell_color_code,
                        name = cell_color
                    )
                } else if (isTRUE(color_as_factor)) {
                    number_colors <- length(unique(factor_data))
                    cell_color_code <- set_default_color_discrete_cell(
                        instrs = instrs
                    )(n = number_colors)
                    names(cell_color_code) <- unique(factor_data)

                    cell_color_code[["other"]] <- other_cell_color
                    pl <- pl +
                        ggplot2::scale_fill_manual(
                            values = cell_color_code,
                            name = cell_color
                        )
                }

                ## plot centers
                if (isTRUE(show_cluster_center) &&
                    (isTRUE(color_as_factor) ||
                        class_cell_color %in% c("character", "factor"))) {
                    pl <- pl + ggplot2::geom_point(
                        data = annotated_DT_centers,
                        aes_string2(
                            x = "center_1",
                            y = "center_2",
                            color = cell_color
                        ),
                        size = center_point_size,
                        shape = 19
                    )
                }

                ## plot labels
                if (isTRUE(show_center_label)) {
                    pl <- pl + ggrepel::geom_text_repel(
                        data = annotated_DT_centers,
                        aes_string2(
                            x = "center_1",
                            y = "center_2",
                            label = cell_color
                        ),
                        size = label_size,
                        fontface = label_fontface
                    )
                }
            }
        }
    }



    ## lastly overlay POINTS ##
    ## first plot other non-selected cells
    if ((!is.null(select_cells) || !is.null(select_cell_groups)) &&
        isTRUE(show_other_cells)) {
        pl <- pl + ggplot2::geom_point(
            data = cell_locations_metadata_other,
            aes_string2(x = sdimx, sdimy),
            color = "black",
            show.legend = FALSE,
            size = other_point_size,
            alpha = point_alpha
        )
    }

    ## plot selected cells
    pl <- pl + ggplot2::geom_point(
        data = cell_locations_metadata_selected,
        aes_string2(x = sdimx, y = sdimy),
        show.legend = FALSE,
        color = "black",
        size = point_size,
        alpha = point_alpha
    )


    return(pl)
}










# spatial in situ ####

#' @title Plot cell polygon layer
#' @name plot_cell_polygon_layer
#' @description Low level function to plot a polygon
#' @return ggplot
#' @details This functions plots a polygon based on spatial cell information.
#' This is most likely a polygon that corresponds to the cell shape.
#' @keywords internal
#' @noRd
plot_cell_polygon_layer <- function(ggobject = NULL,
    instrs = NULL,
    polygon_dt,
    polygon_grouping = "poly_ID",
    fill = NULL,
    poly_fill_gradient = NULL,
    fill_gradient_midpoint = NULL,
    fill_gradient_style = "divergent",
    fill_as_factor = TRUE,
    fill_code = NULL,
    show_legend = TRUE,
    bg_color = "black",
    color = "black",
    alpha = 0.5,
    linewidth = 2) {
    # check fill column
    if (!is.null(fill)) {
        if (isTRUE(fill_as_factor)) {
            polygon_dt[, "final_fill" := as.factor(get(fill))]
        } else {
            polygon_dt[, "final_fill" := get(fill)]
        }
    }

    # create layer
    if (!is.null(ggobject) &&
        methods::is(ggobject, "ggplot")) {
        pl <- ggobject
    } else {
        pl <- ggplot2::ggplot()
    }


    # specific fill color for polygon shapes
    if (!is.null(fill)) {
        pl <- pl + ggplot2::geom_polygon(
            data = polygon_dt,
            aes_string2(
                x = "x",
                y = "y",
                group = polygon_grouping,
                fill = "final_fill"
            ),
            show.legend = show_legend,
            alpha = alpha,
            color = color,
            linewidth = linewidth
        )

        # manual fill colors for factor values
        if (isTRUE(fill_as_factor)) {
            if (!is.null(fill_code)) {
                pl <- pl + ggplot2::scale_fill_manual(values = fill_code)
            } else {
                fill_values_names <- unique(polygon_dt[["final_fill"]])
                fill_code <- set_default_color_discrete_poly(
                    instrs = instrs
                )(length(fill_values_names))
                names(fill_code) <- fill_values_names
                pl <- pl + ggplot2::scale_fill_manual(values = fill_code)
            }
        }

        # gradient fill colors for numerical values
        if (!isTRUE(fill_as_factor)) {
            if (is.null(fill_gradient_midpoint)) {
                fill_gradient_midpoint <-
                    stats::median(polygon_dt[["final_fill"]])
            }

            pl <- pl + set_default_color_continuous_poly(
                colors = poly_fill_gradient,
                instrs = instrs,
                midpoint = fill_gradient_midpoint,
                style = fill_gradient_style,
                guide = ggplot2::guide_colorbar(title = ""),
                type = "fill"
            )
        }
    } else {
        pl <- pl + ggplot2::geom_polygon(
            data = polygon_dt,
            aes_string2(x = "x", y = "y", group = "poly_ID"),
            fill = bg_color,
            alpha = alpha,
            color = color,
            linewidth = linewidth
        )
    }

    return(pl)
}









#' @title plot_feature_points_layer
#' @name plot_feature_points_layer
#' @description low level function to plot a points at the spatial in situ level
#' @return ggplot
#' @details This function can plot multiple features over multiple modalities.
#' These plots can get very big very fast.
#' @keywords internal
#' @noRd
plot_feature_points_layer <- function(ggobject,
    instrs = NULL,
    ext,
    spatial_feat_info,
    feats,
    feats_color_code = NULL,
    feat_shape_code = NULL,
    sdimx = "x",
    sdimy = "y",
    color = "feat_ID",
    shape = "feat",
    point_size = 1.5,
    stroke = NULL,
    show_legend = TRUE,
    plot_method = c("ggplot", "scattermore", "scattermost"),
    expand_counts = FALSE,
    count_info_column = "count",
    jitter = c(0, 0),
    verbose = TRUE) {
    # define plotting method
    plot_method <- match.arg(
        arg = plot_method,
        choices = c("ggplot", "scattermore", "scattermost")
    )

    # data.table vars
    feat_ID <- x <- y <- NULL

    spatial_feat_info_subset <- spatial_feat_info[feat_ID %in% unlist(feats)]

    # expand feature coordinates and/or add jitter to coordinates
    if (isTRUE(expand_counts) | !identical(c(0, 0), jitter)) {
        if (length(jitter) == 1L) jitter <- rep(jitter, 2)
        spatial_feat_info_subset <- expand_feature_info(
            spatial_feat_info = spatial_feat_info_subset,
            expand_counts = expand_counts,
            count_info_column = count_info_column,
            jitter = jitter,
            verbose = verbose
        )
    }

    wrap_msg(
        " --| Plotting ", nrow(spatial_feat_info_subset),
        " feature points"
    )

    if (!is.null(ggobject) & inherits(ggobject, "ggplot")) {
        pl <- ggobject
    } else {
        pl <- ggplot2::ggplot()
    }

    # prepare color vector for scattermost
    if (plot_method == "scattermost") {
        if (!is.null(feats_color_code)) {
            scattermost_color <-
                feats_color_code[spatial_feat_info_subset[["feat_ID"]]]
        } else {
            feats_names <- unique(spatial_feat_info_subset[[color]])
            feats_color_code <- set_default_color_discrete_feat(
                instrs = instrs
            )(length(feats_names))
            names(feats_color_code) <- feats_names
            scattermost_color <-
                feats_color_code[spatial_feat_info_subset[["feat_ID"]]]
        }
    }

    pl <- pl + giotto_point(
        plot_method = plot_method,
        data = spatial_feat_info_subset,
        ext = ext,
        mapping = aes_string2(
            x = sdimx,
            y = sdimy,
            color = color,
            shape = shape
        ),
        size = point_size,
        stroke = stroke,
        show.legend = show_legend,

        # specific for scattermost
        scattermost_xy = spatial_feat_info_subset[, .(x, y)],
        scattermost_color = scattermost_color
    )



    # manually set feature color code
    if (!is.null(feats_color_code)) {
        pl <- pl + ggplot2::scale_color_manual(values = feats_color_code)
    } else {
        feats_names <- unique(spatial_feat_info_subset[[color]])
        feats_color_code <- set_default_color_discrete_feat(
            instrs = instrs
        )(length(feats_names))
        names(feats_color_code) <- feats_names
        pl <- pl + ggplot2::scale_color_manual(values = feats_color_code)
    }

    # manually set feature shape color code
    if (!is.null(feat_shape_code)) {
        pl <- pl + ggplot2::scale_shape_manual(values = feat_shape_code)
    }

    return(pl)
}








#' @title plot_feature_raster_density_layer
#' @name plot_feature_raster_density_layer
#' @description low level function to plot density plots at the spatial
#' in situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
#' @noRd
plot_feature_raster_density_layer <- function(
        ggobject = NULL,
        instrs = NULL,
        spatial_feat_info,
        sel_feat,
        sdimx = "x",
        sdimy = "y",
        alpha = 0.5) {
    # data.table variable
    feat_ID <- NULL

    spatial_feat_info_subset <- spatial_feat_info[feat_ID %in% unlist(sel_feat)]

    if (!is.null(ggobject) & methods::is(ggobject, "ggplot")) {
        pl <- ggobject
    } else {
        pl <- ggplot2::ggplot()
    }

    pl <- pl + ggplot2::stat_density_2d(
        data = spatial_feat_info_subset,
        aes_string2(
            x = sdimx,
            y = sdimy,
            fill = "ggplot2::after_stat(density)"
        ),
        geom = "raster",
        alpha = alpha,
        contour = FALSE
    )
    pl <- pl + labs(title = sel_feat)
    pl <- pl + ggplot2::scale_fill_continuous(type = "viridis")

    return(pl)
}







#' @title plot_feature_hexbin_layer
#' @name plot_feature_hexbin_layer
#' @description low level function to plot hexbins at the spatial in situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
#' @noRd
plot_feature_hexbin_layer <- function(ggobject = NULL,
    instrs = NULL,
    spatial_feat_info,
    sel_feat,
    sdimx = "x",
    sdimy = "y",
    binwidth = NULL,
    min_axis_bins = 10L,
    alpha = 0.5) {
    # data.table variables
    feat_ID <- NULL

    spatial_feat_info_subset <- spatial_feat_info[feat_ID %in% sel_feat]

    # set default binwidth to 1/10 of minor axis
    if (is.null(binwidth)) {
        minorRange <- spatial_feat_info_subset[,
            min(diff(vapply(.SD, range, FUN.VALUE = numeric(2)))),
            .SDcols = c("x", "y")
        ]
        binwidth <- as.integer(minorRange / min_axis_bins)
    }

    if (!is.null(ggobject) & methods::is(ggobject, "ggplot")) {
        pl <- ggobject
    } else {
        pl <- ggplot2::ggplot()
    }

    pl <- pl + ggplot2::geom_hex(
        data = spatial_feat_info_subset,
        aes_string2(x = sdimx, y = sdimy),
        binwidth = binwidth,
        alpha = alpha
    )
    pl <- pl + labs(title = sel_feat)
    return(pl)
}












# image ####

#' @title plot_spat_image_layer_ggplot
#' @name plot_spat_image_layer_ggplot
#' @description create background image in ggplot
#' @param gg_obj ggplot2 object
#' @param gobject giotto object
#' @param gimage a giotto image or a list/vector of giotto images
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name for spatial locations
#' @param polygon_feat_type name for feature type associated with
#' polygon information
#' @param \dots additional args passed to `giotto` `.guess_plot_extent`
#' @return ggplot
#' @keywords internal
#' @noRd
plot_spat_image_layer_ggplot <- function(gg_obj,
    gobject,
    gimage,
    spat_unit = NULL,
    feat_type = NULL,
    spat_loc_name = NULL,
    polygon_feat_type = NULL,
    ...) {
    if (is.null(gobject) || is.null(gimage)) {
        stop("A giotto object and a giotto image need to be provided")
    }

    # prefer extent detection from polygon
    if (!is.null(polygon_feat_type)) spat_unit <- polygon_feat_type

    # allows passing of pre-defined `ext` through ...
    e <- .guess_plot_extent(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type, # only used in detection from points
        spat_loc_name = spat_loc_name,
        prefer = c("polygon", "spatlocs", "points"),
        ...
    )

    # Assign image(s) to plot
    gg_obj <- gg_annotation_raster(ggobj = gg_obj, gimage = gimage, ext = e)

    return(gg_obj)
}









# deconvolution ####

#' @title plot_spat_scatterpie_layer_ggplot
#' @name plot_spat_scatterpie_layer_ggplot
#' @description create scatterpie in ggplot
#' @param ggobject ggplot object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spatial_locations spatial locations
#' @param spatial_enrichment spatial enrichment results
#' @param pie_scale amount to scale the pie size if there is no radius mapping
#' exists
#' @param color color of lines within pie chart
#' @param cell_color_code color code for the cell types
#' @return ggplot
#' @keywords internal
#' @noRd
plot_spat_scatterpie_layer_ggplot <- function(ggobject,
    instrs = NULL,
    sdimx = "sdimx",
    sdimy = "sdimy",
    spatial_locations = NULL,
    spatial_enrichment = NULL,
    pie_scale = 0.35,
    color = NA,
    alpha = 1,
    cell_color_code = NULL) {
    # get cell names
    cell_names <- setdiff(colnames(spatial_enrichment), "cell_ID")

    # combine spatial locations and enrichment results
    combined_spat_enrichm <- data.table::merge.data.table(
        x = spatial_enrichment,
        y = spatial_locations,
        by = "cell_ID"
    )

    # plot scatterpie
    pl <- ggobject
    pl <- pl + scatterpie::geom_scatterpie(
        data = combined_spat_enrichm,
        aes_string2(x = sdimx, y = sdimy),
        pie_scale = pie_scale,
        cols = cell_names,
        color = color,
        alpha = alpha
    )

    ## specificy colors to use
    if (!is.null(cell_color_code)) {
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
    } else {
        number_colors <- length(unique(cell_names))
        cell_color_code <- set_default_color_discrete_cell(
            instrs = instrs
        )(n = number_colors)
        names(cell_color_code) <- unique(cell_names)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
    }

    return(pl)
}





# dim reduction ####

#' @title plot_network_layer_ggplot
#' @name plot_network_layer_ggplot
#' @description Visualize cells in network layer according to dimension
#' reduction coordinates
#' @param ggobject ggplot object
#' @param annotated_network_DT annotated network data.table of selected cells
#' @inheritParams plot_params
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @noRd
plot_network_layer_ggplot <- function(
        ggobject,
        instrs = NULL,
        annotated_network_DT,
        edge_alpha = NULL,
        show_legend = TRUE) {
    from_dims <- grep("from_Dim", colnames(annotated_network_DT), value = TRUE)
    to_dims <- grep("to_Dim", colnames(annotated_network_DT), value = TRUE)



    pl <- ggobject

    if (is.null(edge_alpha)) {
        edge_alpha <- 0.05
        pl <- pl + ggplot2::geom_segment(
            data = annotated_network_DT,
            aes_string2(
                x = from_dims[1],
                y = from_dims[2],
                xend = to_dims[1],
                yend = to_dims[2]
            ),
            alpha = edge_alpha,
            show.legend = show_legend
        )
    } else if (is.numeric(edge_alpha)) {
        pl <- pl + ggplot2::geom_segment(
            data = annotated_network_DT,
            aes_string2(
                x = from_dims[1],
                y = from_dims[2],
                xend = to_dims[1],
                yend = to_dims[2]
            ),
            alpha = edge_alpha,
            show.legend = show_legend
        )
    } else if (is.character(edge_alpha)) {
        if (edge_alpha %in% colnames(annotated_network_DT)) {
            pl <- pl + ggplot2::geom_segment(
                data = annotated_network_DT,
                aes_string2(
                    x = from_dims[1],
                    y = from_dims[2],
                    xend = to_dims[1],
                    yend = to_dims[2],
                    alpha = edge_alpha
                ),
                show.legend = show_legend
            )
        }
    }

    return(pl)
}




#' @title plot_point_layer_ggplot
#' @name plot_point_layer_ggplot
#' @description Visualize cells in point layer according to dimension
#' reduction coordinates
#' @param ggobject ggplot object
#' @inheritParams plot_params
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @noRd
plot_point_layer_ggplot <- function(
        ggobject,
        instrs = NULL,
        annotated_DT_selected,
        annotated_DT_other,
        cell_color = NULL,
        color_as_factor = TRUE,
        cell_color_code = NULL,
        cell_color_gradient = NULL,
        gradient_midpoint = 0,
        gradient_style = "divergent",
        gradient_limits = NULL,
        select_cell_groups = NULL,
        select_cells = NULL,
        point_size = 1,
        point_alpha = 1,
        point_border_col = "black",
        point_border_stroke = 0.1,
        show_cluster_center = FALSE,
        show_center_label = TRUE,
        center_point_size = 4,
        center_point_border_col = "black",
        center_point_border_stroke = 0.1,
        label_size = 4,
        label_fontface = "bold",
        edge_alpha = NULL,
        show_other_cells = TRUE,
        other_cell_color = "lightgrey",
        other_point_size = 0.5,
        show_legend = TRUE) {
    pl <- ggobject



    ## first plot other non-selected cells
    if ((!is.null(select_cells) | !is.null(select_cell_groups)) &
        show_other_cells == TRUE) {
        dims <- grep("Dim.", colnames(annotated_DT_other), value = TRUE)
        pl <- pl + ggplot2::geom_point(
            data = annotated_DT_other,
            aes_string2(x = dims[1], y = dims[2]),
            color = other_cell_color,
            show.legend = FALSE,
            size = other_point_size,
            alpha = point_alpha
        )
    }


    ## order of color
    # 1. if NULL then default to lightblue
    # 2. if character vector
    # 2.1 if length of cell_color is longer than 1 and has colors
    # 2.2 if not part of metadata then suppose its color
    # 2.3 part of metadata
    # 2.3.1 numerical column
    # 2.3.2 factor column or character to factor


    ## point layer
    dims <- grep("Dim.", colnames(annotated_DT_selected), value = TRUE)

    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl + ggplot2::geom_point(
            data = annotated_DT_selected,
            aes_string2(x = dims[1], y = dims[2]),
            color = cell_color,
            show.legend = show_legend,
            size = point_size,
            alpha = point_alpha
        )
        # map color for each cell
    } else if (length(cell_color) > 1) {
        if (is.numeric(cell_color) | is.factor(cell_color)) {
            if (nrow(annotated_DT_selected) != length(cell_color)) {
                stop("\n vector needs to be the same lengths as number
                    of cells \n")
            }
            annotated_DT_selected[["temp_color"]] <- cell_color

            pl <- pl + ggplot2::geom_point(
                data = annotated_DT_selected,
                aes_string2(x = dims[1], y = dims[2], fill = "temp_color"),
                show.legend = show_legend, shape = 21,
                size = point_size,
                color = point_border_col,
                stroke = point_border_stroke,
                alpha = point_alpha
            )
        } else if (is.character(cell_color)) {
            if (!all(.is_color_code(cell_color))) {
                stop("cell_color is not numeric,
                    a factor or vector of colors \n")
            }
            pl <- pl + ggplot2::geom_point(
                data = annotated_DT_selected,
                aes_string2(x = dims[1], y = dims[2]),
                show.legend = show_legend,
                shape = 21,
                fill = cell_color,
                size = point_size,
                color = point_border_col,
                stroke = point_border_stroke,
                alpha = point_alpha
            )
        }
    } else if (is.character(cell_color)) {
        # color by col values
        if (!cell_color %in% colnames(annotated_DT_selected)) {
            if (!.is_color_code(cell_color)) {
                stop(cell_color, " is not a color or a column name \n")
            }
            pl <- pl + ggplot2::geom_point(
                data = annotated_DT_selected,
                aes_string2(x = dims[1], y = dims[2]),
                show.legend = show_legend,
                shape = 21,
                fill = cell_color,
                size = point_size,
                color = point_border_col,
                stroke = point_border_stroke,
                alpha = point_alpha
            )
        } else {
            class_cell_color <- class(annotated_DT_selected[[cell_color]])

            # if len = 1 + color_as_factor = FALSE + data is numeric,
            # assume it means data col index with non-factor data
            if ((class_cell_color == "integer" |
                class_cell_color == "numeric") & color_as_factor == FALSE) {
                # set upper and lower limits
                if (!is.null(gradient_limits) &
                    is.vector(gradient_limits) & length(gradient_limits) == 2) {
                    lower_lim <- gradient_limits[[1]]
                    upper_lim <- gradient_limits[[2]]

                    numeric_data <- annotated_DT_selected[[cell_color]]
                    limit_numeric_data <- ifelse(numeric_data > upper_lim,
                        upper_lim,
                        ifelse(numeric_data < lower_lim,
                            lower_lim, numeric_data
                        )
                    )
                    annotated_DT_selected[[cell_color]] <- limit_numeric_data
                }

                pl <- pl + ggplot2::geom_point(
                    data = annotated_DT_selected,
                    aes_string2(x = dims[1], y = dims[2], fill = cell_color),
                    show.legend = show_legend,
                    shape = 21,
                    size = point_size,
                    color = point_border_col,
                    stroke = point_border_stroke,
                    alpha = point_alpha
                )
            } else {
                # convert character or numeric to factor
                if (color_as_factor == TRUE) {
                    factor_data <- factor(annotated_DT_selected[[cell_color]])
                    annotated_DT_selected[[cell_color]] <- factor_data
                }

                # if you want to show centers or labels then calculate centers
                if (show_cluster_center == TRUE || show_center_label == TRUE) {
                    annotated_DT_centers <- annotated_DT_selected[, .(
                        center_1 = stats::median(get(dims[1])),
                        center_2 = stats::median(get(dims[2]))
                    ), by = cell_color]
                    factor_center_data <- factor(
                        annotated_DT_centers[[cell_color]]
                    )
                    annotated_DT_centers[[cell_color]] <- factor_center_data
                }

                pl <- pl + ggplot2::geom_point(
                    data = annotated_DT_selected,
                    aes_string2(x = dims[1], y = dims[2], fill = cell_color),
                    show.legend = show_legend,
                    shape = 21,
                    size = point_size,
                    color = point_border_col,
                    stroke = point_border_stroke,
                    alpha = point_alpha
                )


                ## plot centers
                if (show_cluster_center == TRUE &
                    (color_as_factor == TRUE |
                        class_cell_color %in% c("character", "factor"))) {
                    pl <- pl + ggplot2::geom_point(
                        data = annotated_DT_centers,
                        aes_string2(
                            x = "center_1",
                            y = "center_2",
                            fill = cell_color
                        ),
                        color = center_point_border_col,
                        stroke = center_point_border_stroke,
                        size = center_point_size,
                        shape = 21,
                        alpha = point_alpha
                    )
                }

                ## plot labels
                if (isTRUE(show_center_label)) {
                    pl <- pl + ggrepel::geom_text_repel(
                        data = annotated_DT_centers,
                        aes_string2(
                            x = "center_1",
                            y = "center_2",
                            label = cell_color
                        ),
                        size = label_size,
                        fontface = label_fontface
                    )
                }
            }


            ## specificy colors to use
            if (!is.null(cell_color_code)) {
                pl <- pl +
                    ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (isTRUE(color_as_factor)) {
                number_colors <- length(unique(factor_data))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instrs
                )(n = number_colors)
                names(cell_color_code) <- unique(factor_data)
                pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
            } else if (!isTRUE(color_as_factor)) {
                if (is.null(gradient_midpoint)) {
                    gradient_midpoint <-
                        stats::median(annotated_DT_selected[[cell_color]])
                }
                pl <- pl + set_default_color_continuous_cell(
                    colors = cell_color_gradient,
                    instrs = instrs,
                    midpoint = gradient_midpoint,
                    style = gradient_style,
                    type = "fill"
                )
            }
        }
    }
    return(pl)
}




#' @title plot_point_layer_ggplot_noFILL
#' @name plot_point_layer_ggplot_noFILL
#' @description Visualize cells in point layer according to dimension
#' reduction coordinates without borders
#' @param ggobject ggplot object
#' @inheritParams plot_params
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @noRd
plot_point_layer_ggplot_noFILL <- function(ggobject,
    plot_method = "ggplot",
    instrs = NULL,
    annotated_DT_selected,
    annotated_DT_other,
    cell_color = NULL,
    color_as_factor = TRUE,
    cell_color_code = NULL,
    cell_color_gradient = NULL,
    gradient_midpoint = 0,
    gradient_style = "divergent",
    gradient_limits = NULL,
    select_cell_groups = NULL,
    select_cells = NULL,
    point_size = 1,
    point_alpha = 1,
    show_cluster_center = FALSE,
    show_center_label = TRUE,
    center_point_size = 4,
    label_size = 4,
    label_fontface = "bold",
    edge_alpha = NULL,
    show_other_cells = TRUE,
    other_cell_color = "lightgrey",
    other_point_size = 0.5,
    show_legend = TRUE,
    ...) {
    pl <- ggobject

    plot_method <- match.arg(plot_method, c("ggplot", "scattermore"))


    ## first plot other non-selected cells
    if ((!is.null(select_cells) |
        !is.null(select_cell_groups)) & show_other_cells == TRUE) {
        dims <- grep("Dim.", colnames(annotated_DT_other), value = TRUE)
        pl <- pl + ggplot2::geom_point(
            data = annotated_DT_other,
            aes_string2(x = dims[1], dims[2]),
            color = other_cell_color,
            show.legend = FALSE,
            size = other_point_size,
            alpha = point_alpha
        )
    }


    ## order of color
    # 1. if NULL then default to lightblue
    # 2. if character vector
    # 2.1 if length of cell_color is longer than 1 and has colors
    # 2.2 if not part of metadata then suppose its color
    # 2.3 part of metadata
    # 2.3.1 numerical column
    # 2.3.2 factor column or character to factor


    ## point layer
    dims <- grep("Dim.", colnames(annotated_DT_selected), value = TRUE)

    if (is.null(cell_color)) {
        cell_color <- "lightblue"
        pl <- pl + giotto_point(
            plot_method = plot_method,
            data = annotated_DT_selected,
            mapping = aes_string2(x = dims[1], dims[2]),
            color = cell_color,
            show.legend = show_legend,
            size = point_size,
            alpha = point_alpha,
            ...
        )
    } else if (length(cell_color) > 1) {
        if (is.numeric(cell_color) | is.factor(cell_color)) {
            if (nrow(annotated_DT_selected) != length(cell_color)) {
                stop("\n vector needs to be the same lengths as
                    number of cells \n")
            }
            annotated_DT_selected[["temp_color"]] <- cell_color

            pl <- pl + giotto_point(
                data = annotated_DT_selected,
                aes_string2(x = dims[1], y = dims[2], color = "temp_color"),
                show.legend = show_legend,
                shape = 19,
                size = point_size,
                alpha = point_alpha
            )
        } else if (is.character(cell_color)) {
            if (!all(.is_color_code(cell_color))) {
                stop("cell_color is not numeric, a factor or vector
                    of colors \n")
            }
            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = annotated_DT_selected,
                mapping = aes_string2(x = dims[1], y = dims[2]),
                show.legend = show_legend,
                shape = 19,
                fill = cell_color,
                size = point_size,
                alpha = point_alpha,
                ...
            )
        }
    } else if (is.character(cell_color)) {
        if (!cell_color %in% colnames(annotated_DT_selected)) {
            if (!.is_color_code(cell_color)) {
                stop(cell_color, " is not a color or a column name \n")
            }
            pl <- pl + giotto_point(
                plot_method = plot_method,
                data = annotated_DT_selected,
                mapping = aes_string2(x = dims[1], y = dims[2]),
                show.legend = show_legend,
                shape = 19,
                color = cell_color,
                size = point_size,
                alpha = point_alpha,
                ...
            )
        } else {
            class_cell_color <- class(annotated_DT_selected[[cell_color]])

            if ((class_cell_color == "integer" |
                class_cell_color == "numeric") & color_as_factor == FALSE) {
                # set upper and lower limits
                if (!is.null(gradient_limits) &
                    is.vector(gradient_limits) & length(gradient_limits) == 2) {
                    lower_lim <- gradient_limits[[1]]
                    upper_lim <- gradient_limits[[2]]

                    numeric_data <- annotated_DT_selected[[cell_color]]
                    limit_numeric_data <- ifelse(numeric_data > upper_lim,
                        upper_lim,
                        ifelse(numeric_data < lower_lim, lower_lim,
                            numeric_data
                        )
                    )
                    annotated_DT_selected[[cell_color]] <- limit_numeric_data
                }

                pl <- pl + giotto_point(
                    plot_method = plot_method,
                    data = annotated_DT_selected,
                    mapping = aes_string2(
                        x = dims[1],
                        y = dims[2],
                        color = cell_color),
                    show.legend = show_legend,
                    shape = 19,
                    size = point_size,
                    alpha = point_alpha,
                    ...
                )
            } else {
                # convert character or numeric to factor
                if (color_as_factor == TRUE) {
                    factor_data <- factor(annotated_DT_selected[[cell_color]])
                    annotated_DT_selected[[cell_color]] <- factor_data
                }

                # if you want to show centers or labels then calculate centers
                if (show_cluster_center == TRUE | show_center_label == TRUE) {
                    annotated_DT_centers <- annotated_DT_selected[, .(
                        center_1 = stats::median(get(dims[1])),
                        center_2 = stats::median(get(dims[2]))
                    ), by = cell_color]
                    factor_center_data <-
                        factor(annotated_DT_centers[[cell_color]])
                    annotated_DT_centers[[cell_color]] <- factor_center_data
                }

                pl <- pl + giotto_point(
                    plot_method = plot_method,
                    data = annotated_DT_selected,
                    mapping = aes_string2(
                        x = dims[1],
                        y = dims[2],
                        color = cell_color),
                    show.legend = show_legend,
                    shape = 19,
                    size = point_size,
                    alpha = point_alpha,
                    ...
                )


                ## plot centers
                if (show_cluster_center == TRUE &
                    (color_as_factor == TRUE |
                        class_cell_color %in% c("character", "factor"))) {
                    pl <- pl + giotto_point(
                        plot_method = plot_method,
                        data = annotated_DT_centers,
                        mapping = aes_string2(
                            x = "center_1",
                            y = "center_2",
                            color = cell_color
                        ),
                        size = center_point_size,
                        shape = 19,
                        alpha = point_alpha,
                        ...
                    )
                }

                ## plot labels
                if (show_center_label == TRUE) {
                    pl <- pl + ggrepel::geom_text_repel(
                        data = annotated_DT_centers,
                        aes_string2(
                            x = "center_1",
                            y = "center_2",
                            label = cell_color
                        ),
                        size = label_size,
                        fontface = label_fontface,
                        alpha = point_alpha
                    )
                }
            }


            ## specificy colors to use
            if (!is.null(cell_color_code)) {
                pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)
            } else if (color_as_factor == TRUE) {
                number_colors <- length(unique(factor_data))
                cell_color_code <- set_default_color_discrete_cell(
                    instrs = instrs
                )(n = number_colors)
                names(cell_color_code) <- unique(factor_data)
                pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)
            } else if (color_as_factor == FALSE) {
                if (is.null(gradient_midpoint)) {
                    gradient_midpoint <- stats::median(
                        annotated_DT_selected[[cell_color]]
                    )
                }
                pl <- pl + set_default_color_continuous_cell(
                    colors = cell_color_gradient,
                    instrs = instrs,
                    midpoint = gradient_midpoint,
                    style = gradient_style,
                    type = "color"
                )
            }
        }
    }
    return(pl)
}






# image ####

# TODO - needs update and completion of overlay option


#' @title addGiottoImageToSpatPlot
#' @name addGiottoImageToSpatPlot
#' @description Add a giotto image to a spatial ggplot object post creation
#' @param spatpl a spatial ggplot object
#' @param gimage a giotto image, see \code{\link{createGiottoImage}}
#' @param layer numeric layer on which to add the giotto image. OR takes 'bg' or
#' 'overlay' as input to designate last (bottom/background) or first
#' (top/overlay)
#' @param alpha (optional) add giotto image to plot with transparency.
#' Numeric. From 0 (transparent) to 1 (fully visible)
#' @returns an updated spatial ggplot object
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#'
#' my_giottoimage <- GiottoClass::createGiottoImage(g,
#'     mg_object = system.file("image.png", package = "GiottoVisuals"),
#'     name = "image"
#' )
#'
#' my_spatplot <- spatPlot2D(g, return_plot = TRUE)
#'
#' addGiottoImageToSpatPlot(spatpl = my_spatplot, gimage = my_giottoimage)
#'
#' @export
addGiottoImageToSpatPlot <- function(spatpl = NULL,
    gimage = NULL,
    layer = c("bg", "overlay"),
    alpha = NULL) {
    layer <- match.arg(arg = layer, choices = c("bg", "overlay"))

    if (is.null(spatpl) | is.null(gimage)) {
        stop("A spatial ggplot object and a giotto image need to be given")
    }

    # extract min and max from object
    my_xmax <- gimage@minmax[1]
    my_xmin <- gimage@minmax[2]
    my_ymax <- gimage@minmax[3]
    my_ymin <- gimage@minmax[4]

    # convert giotto image object into array
    img_array <- as.numeric(gimage@mg_object[[1]])

    # add transparency if needed
    if (!is.null(alpha) & is.numeric(alpha)) {
        img_array <- add_img_array_alpha(
            x = img_array,
            alpha = alpha
        )
    }

    # extract adjustments from object
    xmax_b <- gimage@boundaries[1]
    xmin_b <- gimage@boundaries[2]
    ymax_b <- gimage@boundaries[3]
    ymin_b <- gimage@boundaries[4]

    newpl <- spatpl + annotation_raster(
        img_array,
        xmin = my_xmin - xmin_b,
        xmax = my_xmax + xmax_b,
        ymin = my_ymin - ymin_b,
        ymax = my_ymax + ymax_b
    )

    # position new layer
    if (layer == "bg") {
        # move image to background
        nr_layers <- length(newpl$layers)
        newpl$layers <- c(
            newpl$layers[[nr_layers]],
            newpl$layers[seq_len(nr_layers - 1)]
        )
    } else if (layer == "overlay") {} # keep image on top


    return(newpl)
}
