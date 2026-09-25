#' @name gg_annotation_raster
#' @title Append image to ggplot as annotation_raster
#' @description
#' Add a spatially mapped image to a *ggplot2* `gg` object.
#' For terra-based images, the image will be a cropped and sampled version
#' of the full size image on disk that has sufficient resolution for the size
#' of the plot requested.
#'
#' @param ggobj ggplot2 `gg` object
#' @param gimage `giottoLargeImage`, `giottoImage` or `list` thereof
#' @param ext Object that responds to `ext()`. Defines the plot spatial ROI
#' This extent defines which portions of the image(s) will be plotted/should
#' be sampled for. The default is the same extent as the image.
#' @param geom_blank logical. Whether to apply `[ggplot2::geom_blank()]` to the
#' `gg` object so that the image can be plotted by itself.
#' @param \dots additional params to pass
#' @details
#' No ... params are implemented for `giottoImage`. \cr ... params for
#' `giottoLargeImage` and `giottoAffineImage` pass to `?auto_image_resample`
#' @return `gg` object with images to plot appended as annotation rasters
#' @examples
#' gimg <- GiottoData::loadSubObjectMini("giottoLargeImage")
#' gg <- ggplot2::ggplot()
#' out <- GiottoVisuals::gg_annotation_raster(gg, gimg)
#' print(out)
NULL

# * list ####
#' @rdname gg_annotation_raster
#' @export
setMethod(
    "gg_annotation_raster",
    signature(ggobj = "ggUnionClass", gimage = "list"),
    function(ggobj, gimage, ext = NULL, geom_blank = TRUE, ...) {
        # apply geom_blank
        ext <- ext %null% ext(gimage[[1L]])
        if (geom_blank) ggobj <- .gg_geom_blank(ggobj, ext)

        # attach images in a loop
        for (i in seq_along(gimage)) {
            ggobj <- gg_annotation_raster(
                ggobj, gimage[[i]],
                ext = ext,
                geom_blank = FALSE, # hardcode FALSE since already done.
                ...
            )
        }
        return(ggobj)
    }
)

# * giottoImage ####
#' @rdname gg_annotation_raster
#' @export
setMethod(
    "gg_annotation_raster",
    signature(ggobj = "ggUnionClass", gimage = "giottoImage"),
    function(ggobj, gimage, ext = NULL, geom_blank = TRUE, ...) {
        # apply geom_blank
        ext <- ext %null% ext(gimage)
        if (geom_blank) ggobj <- .gg_geom_blank(ggobj, ext)

        # convert giotto image object into array
        img_array <- as.numeric(gimage@mg_object[[1]])

        # append to ggobj
        ggobj <- .gg_append_imagearray(ggobj, img_array, ext)

        # TODO geom_raster to accommodate single-channel
        return(ggobj)
    }
)

# * giottoLargeImage ####
#' @rdname gg_annotation_raster
#' @export
setMethod(
    "gg_annotation_raster",
    signature(ggobj = "ggUnionClass", gimage = "giottoLargeImage"),
    function(ggobj, gimage, ext = NULL, geom_blank = TRUE, ...) {
        # geom_blank
        ext <- ext %null% ext(gimage)
        if (geom_blank) ggobj <- .gg_geom_blank(ggobj, ext)

        # resample from extent
        gimage <- .auto_resample_gimage(
            img = gimage,
            plot_ext = ext,
            sample_fun = .sample_gimage,
            ...
        )

        # append raster to gg
        ggobj <- .gg_append_spatraster(ggobj = ggobj, gimage = gimage)

        return(ggobj)
    }
)

# * giottoAffineImage ####
#' @rdname gg_annotation_raster
#' @export
setMethod(
    "gg_annotation_raster",
    signature(ggobj = "ggUnionClass", gimage = "giottoAffineImage"),
    function(ggobj, gimage, ext = NULL, geom_blank = TRUE, ...) {
        # geom_blank
        ext <- ext %null% ext(gimage)
        if (geom_blank) ggobj <- .gg_geom_blank(ggobj, ext)

        # resample from extent
        gimage <- .auto_resample_gimage(
            img = gimage,
            plot_ext = ext,
            sample_fun = .sample_gaffimage,
            ...
        )

        # append raster to gg
        ggobj <- .gg_append_spatraster(ggobj = ggobj, gimage = gimage)

        return(ggobj)
    }
)




# Internals ####

# returns the spatial extent needed for the plot
# ... passes to ext() `giotto` method
.guess_plot_extent <- function(
        gobject, spat_unit = NULL, spat_loc_name = NULL, ext = NULL, ...) {
    if (!is.null(ext)) ext <- ext(ext) # normalize to `SpatExtent` class
    # if ext already given, directly return
    if (inherits(ext, "SpatExtent")) {
        return(ext)
    }

    # find extent from one of poly, spatlocs, points, in that order of pref
    e <- ext(
        gobject,
        spat_unit = spat_unit,
        all_data = FALSE,
        verbose = FALSE,
        name = list(spatlocs = spat_loc_name),
        # `name` only passes to `getSpatialLocations()` if spatlocs are
        # present and used to find extent
        ... # You can ensure they are used by setting prefer = "spatlocs
    )

    if (is.null(e)) {
        stop(wrap_txt(
            "No `ext` provided and no spatial locations or polygons discovered.
            Cannot determine largeImage resample extent"
        ))
    }
    return(e)
}

# internal to convert a SpatExtent into a data.frame with x and y values that
# ggplot2 can use to determine bounds of placement
.ext_to_dummy_df <- function(x) {
    data.frame(
        sdimx = x[][c(1, 2)],
        sdimy = x[][c(3, 4)],
        row.names = NULL
    )
}

# apply a region to plot to the gg object. Input should be a SpatExtent or
# coercible. Returns ggobject with geom_blank assigned
.gg_geom_blank <- function(ggobj, e) {
    # NSE vars
    sdimx <- sdimy <- NULL

    # create minimal dummy value data.frame of spatial locations that cover
    # the spatial region to plot
    bounds_dt <- .ext_to_dummy_df(e)
    # assign region to plot
    ggobj <- ggobj + geom_blank(data = bounds_dt, aes(sdimx, sdimy))
    return(ggobj)
}

#' @name auto_image_resample
#' @title Optimized image resampling
#' @description
#' Downsample terra-based images for plotting. Uses \code{\link[terra]{window}}
#' to set a virtual reading window on the image before calling
#' \code{\link[terra]{spatSample}}, so only the relevant spatial ROI is read
#' from disk. This avoids materializing a crop to disk for large ROIs while
#' still focusing sampling on the plot region.
#' @param img giotto image to plot
#' @param plot_ext extent of plot (defaults to the image extent)
#' @param img_border numeric. Default = 0.125. If greater than 0, expand
#' `plot_ext` by this fraction on each side before setting the window. See
#' details.
#' @param max_sample numeric. Default = 5e5. Maximum number of values to sample
#' from the image. Globally settable with option "giotto.plot_img_max_sample"
#' @details
#' **img_border**
#' Expands the window extent used for reading the image. This prevents the
#' image from being cut off sharply at the plot boundary, since plot extents
#' are typically defined by centroids and polygons may hang over the edge.
#' @returns a giotto image resampled within the plot window
#' @examples
#' \dontrun{
#' img <- GiottoData::loadSubObjectMini("giottoLargeImage")
#' .auto_resample_gimage(img)
#' }
#' @seealso \code{\link[terra]{window}}, \code{\link[terra]{spatSample}}
#' @keywords internal
.auto_resample_gimage <- function(img,
    plot_ext = NULL,
    img_border = 0.125,
    sample_fun = .sample_gimage,
    max_sample = getOption("giotto.plot_img_max_sample", 5e5)) {
    # determine ext to use
    crop_ext <- if (is.null(plot_ext)) {
        ext(img) # fallback to img extent
    } else {
        plot_ext <- ext(plot_ext) # if ext specified
        if (img_border > 0) { # apply border expansion
            plot_ext <- plot_ext |>
                as.polygons() |>
                rescale(1 + img_border) |>
                ext()
        }
        # normalize extent when larger than available
        terra::intersect(ext(img), plot_ext) # NULL if no intersect
    }

    # check image is in plot_ext
    if (is.null(crop_ext)) {
        warning(sprintf("image '%s' is not within the plotting window", 
            objName(img)), call. = FALSE)
        return(NULL)
    }

    img <- crop(img, crop_ext)
    sample_fun(img, size = max_sample)
}

# Pull sampled values from original image into target spatial mapping
# Returns a `giottoLargeImage`
.sample_gimage <- function(x, size) {
    # "display" lets GDAL read from overviews when present. "regular" ignores
    # them and costs ~2s per call on a JPEG-2000 VRT.
    x[] <- terra::spatSample(x[],
        size = size,
        method = "display"
    )
    return(x)
}

.sample_gaffimage <- function(x, size) {
    res <- x@funs$realize_magick(size = size)
    return(res)
}




# make an image array compatible with ggplot::annotation_raster()
# maxval is the cutoff after which everything is max intensity
# returns: raster
.gg_imgarray_2_raster <- function(x, maxval = NULL, col = NULL) {
    nlyr <- dim(x)[3L] # number of channels/layers
    if (is.na(nlyr)) nlyr <- 1L
    # NOTE: 4 layers allowed (rgba), but may conflict with actual 4 info
    # layer cases which SHOULD be converted to 3 layer
    #
    # more than 4 layers -> directly ignore layers past the 3rd
    if (nlyr > 4L) {
        nlyr <- 3L
        x <- x[, , seq_len(3)]
    }

    # handle NaN values -- set as max value of that layer
    # these may arise due to save artefacting when values are larger than
    # expected
    for (lyr in seq_len(nlyr)) {
        if (is.nan(max(x[, , lyr]))) {
            x[, , lyr][is.nan(x[, , lyr])] <-
                max(x[, , lyr], na.rm = TRUE)
        }
    }

    # handle NA values -- set as 0
    x[is.na(x)] <- 0

    if (nlyr == 1L) {
        # SINGLE CHANNEL #
        # max window cutoff
        maxval <- maxval %na% quantile(x, 0.99)
        if (!is.null(maxval)) x[x > maxval] <- maxval
        # colorize
        if (is.null(col)) {
            col <- getMonochromeColors("white", n = 256)
        }
        r <- .colorize_single_channel_raster(x, col = col)
    } else {
        # RGB EXPECTED #
        # convert to range 0:1 (needed for as.raster())
        x <- scales::rescale(x, to = c(0, 1))
        r <- grDevices::as.raster(x)
    }

    return(r)
}




# `x` is array to use
# `col` is character vector of colors to use
.colorize_single_channel_raster <- function(x, col) {
    if (!is.na(dim(x)[3L])) x <- x[, , 1L] # convert to matrix
    r <- range(x, na.rm = TRUE)
    x <- (x - r[1]) / (r[2] - r[1])
    x <- round(x * (length(col) - 1) + 1)
    x[] <- col[x]
    terra::as.raster(x)
}

# append image array to a gg object
.gg_append_imagearray <- function(ggobj, a, ext) {
    # append to ggobj
    extent <- ext(ext)[seq_len(4L)]
    ggobj <- ggobj + annotation_raster(a,
        xmin = extent[["xmin"]], xmax = extent[["xmax"]],
        ymin = extent[["ymin"]], ymax = extent[["ymax"]]
    )
}

# append a giotto image object containing a SpatRaster that has already been
# resampled/pulled into memory. Output is a `gg` object
.gg_append_spatraster <- function(ggobj, gimage) {
    if (is.null(gimage)) return(ggobj) # passthrough

    # convert gimage to a raster
    a <- terra::as.array(gimage@raster_object) %>%
        .gg_imgarray_2_raster(
            maxval = gimage@max_window,
            col = gimage@colors
        )

    ggobj <- .gg_append_imagearray(ggobj, a, ext(gimage))
    return(ggobj)
}
