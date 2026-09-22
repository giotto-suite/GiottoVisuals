# A spatial network is stored as an igraph carrying only vertex names, so the
# endpoint coordinates every drawing path needs are looked up in the spatial
# locations at draw time. These pin down both halves of that: that the lookup
# happens at all, and that an edge with no drawable endpoint is dropped.
#
# The first test is the one that would have caught the regression these were
# written for. A plot with `show_network = TRUE` CONSTRUCTS fine without the
# coordinate columns -- ggplot resolves aesthetics lazily -- and only fails when
# something renders it. Build the plot, do not just return it.

skip_if_no_mini <- function() skip_if_not_installed("GiottoData")

.net_mini <- function() {
    suppressMessages(GiottoData::loadGiottoMini("vizgen", verbose = FALSE))
}

.net_plot <- function(g, ...) {
    suppressMessages(spatPlot2D(g,
        spat_unit = "aggregate", cell_color = "leiden_clus",
        show_network = TRUE, spatial_network_name = "Delaunay_network",
        show_plot = FALSE, save_plot = FALSE, return_plot = TRUE, ...))
}

# layer 1 is the network segments, layer 2 the cells
.edges <- function(p) p$layers[[1L]]$data
.cells <- function(p) p$layers[[2L]]$data


test_that("a network plot renders, not merely builds", {
    skip_if_no_mini()
    p <- .net_plot(.net_mini())
    expect_named(.edges(p),
        c("to", "from", "distance", "weight", "sdimx_begin", "sdimy_begin",
            "sdimx_end", "sdimy_end"),
        ignore.order = TRUE)
    expect_no_error(ggplot2::ggplot_build(p))
})

test_that("a view narrows the edges with the cells", {
    skip_if_no_mini()
    g <- .net_mini()
    e <- GiottoClass::ext(g, spat_unit = "aggregate", prefer = "polygon")
    g <- GiottoClass::crop(g,
        terra::ext(c(e[1], mean(e[1:2]), e[3], mean(e[3:4]))), view = "roi")

    full <- .net_plot(g)
    roi <- .net_plot(g, view = "roi")

    expect_lt(nrow(.cells(roi)), nrow(.cells(full)))
    expect_lt(nrow(.edges(roi)), nrow(.edges(full)))
    # the property that matters: nothing is drawn to a cell that is not
    expect_length(
        setdiff(unique(c(.edges(roi)$from, .edges(roi)$to)),
            .cells(roi)$cell_ID),
        0L)
    expect_no_error(ggplot2::ggplot_build(roi))
})

test_that("a space moves the edges with the locations", {
    skip_if_no_mini()
    # Nothing in the network path is space-aware; this passes because the
    # locations it joins to were already transformed.
    g <- GiottoClass::spatShift(.net_mini(), dx = 1000, space = "shifted")
    a <- .edges(.net_plot(g))
    b <- .edges(.net_plot(g, space = "shifted"))

    expect_identical(nrow(a), nrow(b))
    expect_equal(b$sdimx_begin - a$sdimx_begin, rep(1000, nrow(a)))
    expect_equal(b$sdimx_end - a$sdimx_end, rep(1000, nrow(a)))
    expect_equal(a$sdimy_begin, b$sdimy_begin)
})

# The join rule itself (inner, and what it does to an absent endpoint) is
# GiottoClass's; see its test-spatial-network-annotate.R.


# window / view composition ####################################################

test_that("an in-situ window composes with the view rather than replacing it", {
    skip_if_no_mini()
    g <- suppressMessages(GiottoData::loadGiottoMini("vizgen", verbose = FALSE))
    e <- GiottoClass::ext(g, spat_unit = "aggregate", prefer = "polygon")
    # view keeps the left half, window keeps the bottom half
    g <- GiottoClass::crop(g,
        terra::ext(c(e[1], mean(e[1:2]), e[3], e[4])), view = "left")
    win <- c(e[3], mean(e[3:4]))

    n <- function(...) {
        p <- suppressMessages(spatInSituPlotPoints(g,
            polygon_feat_type = "aggregate", polygon_fill = "nr_feats",
            polygon_fill_as_factor = FALSE, show_plot = FALSE,
            save_plot = FALSE, return_plot = TRUE, ...))
        nrow(p$layers[[1L]]$data)
    }

    view_only <- n(view = "left")
    window_only <- n(ylim = win)
    both <- n(view = "left", ylim = win)

    # an intersection: strictly smaller than either on its own. If the window
    # replaced the view (or vice versa) one of these would be an equality.
    expect_lt(both, view_only)
    expect_lt(both, window_only)
})
