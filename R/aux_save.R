#' @title Plot saving
#' @name plot_save
#' @description
#' Functions to save plots to directory of interest.
#' `all_plots_save_function()` is used for plot saving operations. `gpsparam()`
#' is used to generate a set of save parameters and filepath based on available
#' parameter and `giottoInstructions` values.
#' @param gobject giotto object or giottoInstructions
#' @param plot_object ggplot object to plot
#' @param save_dir directory to save to
#' @param save_folder folder in save_dir to save to
#' @param save_name name of plot
#' @param default_save_name default name to save a plot
#' @param save_format format (e.g. png, tiff, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @param ncol number of columns for multiplots
#' @param nrow number of rows for multiplot
#' @param scale scale of plots
#' @param base_width width of plot
#' @param base_height height of plot
#' @param base_aspect_ratio aspect ratio of plot
#' @param units plotting units (e.g. in)
#' @param dpi dpi for each plot if plot is in raster format
#' @param limitsize When TRUE (the default), ggsave will not save images larger
#' than 50x50 inches, to prevent the common error of specifying dimensions in
#' pixels.
#' @param plot_count manually set the plot count that is appended to a
#' default_save_name
#' @param GPSPARAM `giotto_plot_save_param` object. If provided, will be
#' used instead of most other general save params. (save_dir, save_folder,
#' save_name, default_save_name, save_format, base_width, base_height,
#' base_aspect_ratio, units, dpi, plot_count)
#' @param \dots additional parameters to pass downstream save functions.
#' [cowplot::save_plot()] is used for `ggplot2` plots. grDevices png, tiff
#' svg, pdf is used for base and general saving
#' @returns `all_plots_save_function` returns a plot file. `gpsparam` returns
#' a `giotto_plot_save_param` object
#' @seealso \code{\link{showSaveParameters}} \code{\link[cowplot]{save_plot}}
#' \code{\link[grDevices]{png}}
#' \code{\link[grDevices]{tiff}}
#' \code{\link[grDevices]{pdf}}
#' \code{\link[grDevices]{svg}}
NULL





#' @rdname plot_save
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#' df <- data.frame(x = rnorm(5), y = rnorm(5))
#' g_plot <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
#'     ggplot2::geom_point()
#' all_plots_save_function(g, g_plot)
#'
#' @export
all_plots_save_function <- function(
        gobject,
        plot_object,
        save_dir = NULL,
        save_folder = NULL,
        save_name = NULL,
        default_save_name = "giotto_plot",
        save_format = NULL,
        show_saved_plot = FALSE,
        ncol = 1,
        nrow = 1,
        scale = 1,
        base_width = NULL,
        base_height = NULL,
        base_aspect_ratio = NULL,
        units = NULL,
        dpi = NULL,
        limitsize = TRUE,
        plot_count = NULL,
        GPSPARAM = NULL,
        ...) {
    # get save params
    if (is.null(GPSPARAM)) {
        type <- "general"
        if (any("ggplot" %in% class(plot_object))) type <- "gg"
        if (any("plotly" %in% class(plot_object))) type <- "plotly"

        a <- .grab_gpsparam_args()
        GPSPARAM <- do.call(gpsparam, args = c(
            a, list(instructions = instructions(gobject), type = type)
        ))
    }

    checkmate::assert_class(GPSPARAM, "giotto_plot_save_param")

    if (identical(getOption("giotto.verbose"), "debug")) {
        print(GPSPARAM)
    }

    # perform save
    if (any("ggplot" %in% class(plot_object))) {
        .ggplot_save_function(
            gobject = gobject,
            plot_object = plot_object,
            show_saved_plot = show_saved_plot,
            ncol = ncol,
            nrow = nrow,
            scale = scale,
            limitsize = limitsize,
            GPSPARAM = GPSPARAM,
            ...
        )
    } else {
        .general_save_function(
            gobject = gobject,
            plot_object = plot_object,
            show_saved_plot = show_saved_plot,
            GPSPARAM = GPSPARAM,
            ...
        )
    }
}



#' @rdname plot_save
#' @param instructions `giotto` or `giottoInstructions` object
#' @param type `character`. One of `"gg"', '"plotly"', '"general"` to designate
#'   which type of plot to save. This affects which types of outputs are
#'   possible.
#' @export
gpsparam <- function(instructions,
    type = c("gg", "plotly", "general"),
    save_dir = NULL,
    save_folder = NULL,
    save_name = NULL,
    default_save_name = "giotto_plot",
    save_format = NULL,
    dpi = NULL,
    base_width = NULL,
    base_height = NULL,
    base_aspect_ratio = NULL,
    units = NULL,
    plot_count = NULL,
    ... # ignored
) {
    if (!inherits(instructions, c("giotto", "giottoInstructions"))) {
        stop(
            "`instructions` must be either a `giotto` or",
            "`giottoInstructions` object."
        )
    }
    instrs <- instructions # shortname
    checkmate::assert_character(type)
    if (!length(type) == 1L) {
        stop("Single `type` must be specified.")
    }

    ## save format -------------------------------------------------------- ##
    save_format <- save_format %null%
        instructions(instrs, param = "plot_format")

    save_format <- switch(type,
        "gg" = save_format,
        "plotly" = "html",
        "general" = match.arg(save_format, c("png", "tiff", "pdf", "svg"))
    )

    ## get save information and set defaults ------------------------------ ##
    save_dir <- save_dir %null% instructions(instrs, param = "save_dir")
    custom_plot_count <- is.null(plot_count)
    plot_count <- plot_count %null% getOption("giotto.plot_count", 1)
    dpi <- dpi %null% instructions(instrs, param = "dpi")
    base_width <- base_width %null% instructions(instrs, param = "width")
    base_height <- base_height %null% instructions(instrs, param = "height")
    base_aspect_ratio <- base_aspect_ratio %null% 1.1
    units <- units %null% instructions(instrs, param = "units")


    ## checking ----------------------------------------------------------- ##
    dpi <- as.numeric(dpi)
    base_width <- as.numeric(base_width)
    base_height <- as.numeric(base_height)
    base_aspect_ratio <- as.numeric(base_aspect_ratio)
    if (is.na(save_dir)) save_dir <- getwd()


    # build filepath ------------------------------------------------------ ##
    if (is.null(save_name)) {
        save_name <- default_save_name
        save_name <- paste0(plot_count, "-", save_name)
        if (custom_plot_count) {
            on.exit(options("giotto.plot_count" = plot_count + 1L), # increment
                add = TRUE
            )
        }
    }

    if (!is.null(save_folder)) {
        file_location <- file.path(save_dir, save_folder)
    } else {
        file_location <- save_dir
    }

    filename <- paste0(save_name, ".", save_format)
    fullpath <- file.path(file_location, filename)

    # create params object ------------------------------------------------ ##
    structure(
        list(
            fullpath = fullpath,
            save_format = save_format,
            dpi = dpi,
            base_width = base_width,
            base_height = base_height,
            base_aspect_ratio = base_aspect_ratio,
            units = units
        ),
        class = "giotto_plot_save_param"
    )
}




#' @title showSaveParameters
#' @name showSaveParameters
#' @description Description of Giotto saving options,
#' links to \code{\link{all_plots_save_function}}
#' @returns Instructions on how to use the automatic plot saving options
#' within Giotto
#' @export
#' @examples
#' showSaveParameters()
showSaveParameters <- function() {
    message("This is a simple guide to help you automatically saving plots. \n")
    message("Importantly, defaults for all these parameters can be set at
    the beginning with createGiottoInstructions() \n")
    message("See
    https://rubd.github.io/Giotto/articles/instructions_and_plotting.html
    for more information and examples \n \n")

    message("Each plotting function in Giotto has 4 important parameters for
    showing and/or saving a plot: \n
        - show_plot: TRUE or FALSE, show the plot to the console
        - return_plot: TRUE or FALSE, return the plot to the console
        (e.g. to further modify or save the plot
        - save_plot: TRUE or FALSE, automatically save the plot
        - save_param: a list of parameters that can be set \n")

    message("The following list of parameters can be provided to save_param: \n
        - save_dir: directory to save the plot to
        - save_folder: if not NULL, a subfolder within save_dir that will
        be created to save the plot to
        - save_name: name of the plot (no extension needed, see save_format)
        - save_format: picture format to use, default is .png
        - ncol: number of columns for multiplots
        - nrow: number of rows for multiplot
        - scale: scale of plots
        - base_width: width of plot
        - base_height: height of plot
        - base_aspect_ratio: ratio of plot
        - units: plotting units (e.g. in)
        - dpi: dpi for each plot if plot is in raster format\n")

    message("Example: \n
        plotfunction(...,
                    save_plot = TRUE,
                    save_param = list(save_name = 'favorite_name',
                    units = 'png'))")
}




# internals ####

# GPSPARAM should be a `giotto_plot_save_param` object if provided
#' @noMd
#' @keywords internal
.ggplot_save_function <- function(gobject,
    plot_object,
    show_saved_plot = FALSE,
    ncol = 1,
    nrow = 1,
    scale = 1,
    limitsize = TRUE,
    GPSPARAM = NULL,
    ...) {
    if (is.null(plot_object)) {
        stop("\t there is no object to plot \t")
    }

    sparam <- GPSPARAM

    # create saving location
    fullpath <- sparam$fullpath
    filename <- basename(fullpath)
    path <- dirname(fullpath)
    save_format <- sparam$save_format

    if (!file.exists(path)) {
        dir.create(path, recursive = TRUE)
    }


    cowplot::save_plot(
        plot = plot_object,
        filename = filename,
        path = path,
        ncol = ncol,
        nrow = nrow,
        scale = scale,
        device = save_format,
        limitsize = limitsize,
        # save param items
        dpi = sparam$dpi,
        units = sparam$units,
        base_width = sparam$base_width,
        base_height = sparam$base_height,
        base_aspect_ratio = sparam$base_aspect_ratio,
        ...
    )

    # show saved plot if requested
    if (isTRUE(show_saved_plot)) {
        if (save_format == "png") {
            if (package_check("png", on_missing = "message")) {
                img <- png::readPNG(source = fullpath)
                grid::grid.raster(img)
            }
        } else if (save_format == "tiff") {
            if (package_check("tiff", on_missing = "message")) {
                img <- tiff::readTIFF(source = fullpath)
                grid::grid.raster(img)
            }
        } else {
            warning("\t only png & tiff are currently supported \t")
        }
    }
}



# GPSPARAM should be a `giotto_plot_save_param` object if provided
#' @noMd
#' @keywords internal
.general_save_function <- function(gobject,
    plot_object,
    show_saved_plot = FALSE,
    GPSPARAM = NULL,
    ...) {
    if (is.null(plot_object)) {
        stop("\t there is no object to plot \t")
    }

    sparam <- GPSPARAM

    fullpath <- sparam$fullpath
    save_format <- sparam$save_format
    dpi <- sparam$dpi
    units <- sparam$units
    base_width <- sparam$base_width
    base_height <- sparam$base_height


    # create saving location
    path <- dirname(fullpath)
    if (!file.exists(path)) {
        dir.create(path, recursive = TRUE)
    }

    if (any("plotly" %in% class(plot_object))) {
        htmlwidgets::saveWidget(
            plotly::as_widget(plot_object),
            file = fullpath
        )
    } else {
        switch(save_format,
            "png" = {
                grDevices::png(
                    filename = fullpath, width = base_width,
                    height = base_height, res = dpi, units = units, ...
                )
                print(plot_object)
                grDevices::dev.off()
            },
            "tiff" = {
                grDevices::tiff(
                    filename = fullpath, width = base_width,
                    height = base_height, units = units, ...
                )
                print(plot_object)
                grDevices::dev.off()
            },
            "pdf" = {
                grDevices::pdf(
                    file = fullpath, width = base_width,
                    height = base_height, useDingbats = FALSE, ...
                )
                print(plot_object)
                grDevices::dev.off()
            },
            "svg" = {
                grDevices::svg(
                    filename = fullpath, width = base_width,
                    height = base_height, ...
                )
                print(plot_object)
                grDevices::dev.off()
            }
        )

        # show saved plot if requested
        if (isTRUE(show_saved_plot)) {
            switch(save_format,
                "png" = {
                    if (package_check("png", on_missing = "message")) {
                        img <- png::readPNG(source = fullpath)
                        grid::grid.raster(img)
                    }
                },
                "tiff" = {
                    if (package_check("tiff", on_missing = "message")) {
                        img <- tiff::readTIFF(source = fullpath)
                        grid::grid.raster(img)
                    }
                },
                warning("\t only png & tiff are currently supported \t")
            )
        }
    }
}


# get expected save params from one stack frame up.
.grab_gpsparam_args <- function() {
    expected_save_argnames <- c(
        "save_dir", "save_folder", "save_name", "default_save_name",
        "save_format", "dpi", "base_width", "base_height", "base_aspect_ratio",
        "units", "plot_count"
    )

    get_args_list(toplevel = 2L, keep = expected_save_argnames)
}


#' @export
print.giotto_plot_save_param <- function(x, ...) {
    cat(sprintf("<%s>\n", class(x)))
    print_list(x)
}

# GPSPARAM should be a `giotto_plot_save_param`
.plot_px_area <- function(GPSPARAM) {
    dims <- c(GPSPARAM$base_height, GPSPARAM$base_width)
    pxdims <- switch(GPSPARAM$units,
        "in" = dims * GPSPARAM$dpi,
        "cm" = (dims / 2.54) * GPSPARAM$dpi,
        "mm" = (dims / 25.4) * GPSPARAM$dpi,
        "px" = dims
    )
    round(prod(pxdims))
}
