# plotLabelProportions(): stacked bars of a label column within a group
# column, from calculateLabelProportions(), a stored composition enrichment,
# or a precomputed table. Ground truth is the counts the fixture is built
# from.

.lp_gobject <- function(n_cells = 120L, seed = 3L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(10L, n_cells, density = 0.5,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(10L))
    colnames(m) <- paste0("c", seq_len(n_cells))
    locs <- data.frame(sdimx = stats::runif(n_cells, 0, 100),
        sdimy = stats::runif(n_cells, 0, 100), cell_ID = colnames(m))
    g <- GiottoClass::createGiottoObject(expression = m,
        spatial_locs = locs, verbose = FALSE)
    # group a: 30 x, 10 y; group b: 10 x, 30 y; group c: 20 x, 20 y
    grp <- rep(c("a", "b", "c"), each = 40L)
    lab <- c(rep(c("x", "y"), c(30L, 10L)), rep(c("x", "y"), c(10L, 30L)),
        rep(c("x", "y"), c(20L, 20L)))
    GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(cell_ID = colnames(m), grp = grp,
            lab = lab))
}

.lp <- function(...) {
    plotLabelProportions(..., show_dendrogram = FALSE,
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
}

.lp_value <- function(p, g, l) {
    p$data$value[p$data$group == g & p$data$label == l]
}

test_that("labels give the share of each group's cells", {
    g <- .lp_gobject()
    p <- .lp(g, groups = "grp", labels = "lab")
    expect_s3_class(p, "ggplot")
    expect_equal(.lp_value(p, "a", "x"), 0.75)
    expect_equal(.lp_value(p, "b", "y"), 0.75)
    expect_equal(.lp_value(p, "c", "x"), 0.5)
    sums <- tapply(p$data$value, p$data$group, sum)
    expect_equal(as.vector(sums), rep(1, 3))
})

test_that("value = 'count' gives cells", {
    g <- .lp_gobject()
    p <- .lp(g, groups = "grp", labels = "lab", value = "count")
    expect_equal(.lp_value(p, "a", "x"), 30)
    expect_equal(as.vector(tapply(p$data$value, p$data$group, sum)),
        rep(40, 3))
})

test_that("a precomputed table draws the same bars", {
    g <- .lp_gobject()
    tab <- GiottoClass::calculateLabelProportions(g, labels = "lab",
        group_method = "table", groups = "grp", output = "data.table",
        verbose = FALSE)
    p1 <- .lp(g, groups = "grp", labels = "lab")
    p2 <- .lp(g, proportions = tab)
    expect_equal(p2$data[order(group, label), value],
        p1$data[order(group, label), value])
    expect_identical(p2$labels$y, "grp")
})

test_that("an enrichment is averaged within each group", {
    g <- .lp_gobject()
    g <- GiottoClass::createSpatialNetwork(g, method = "kNN", k = 5L,
        name = "knn", verbose = FALSE)
    g <- GiottoClass::calculateLabelProportions(g, labels = "lab",
        group_method = "spatialnetwork", spatial_network_name = "knn",
        name = "nb", output = "gobject", verbose = FALSE)
    p <- .lp(g, groups = "grp", spat_enr_name = "nb")
    enr <- GiottoClass::getSpatialEnrichment(g, name = "nb",
        output = "data.table")
    meta <- GiottoClass::pDataDT(g)
    a_ids <- meta$cell_ID[meta$grp == "a"]
    expect_equal(.lp_value(p, "a", "x"),
        mean(enr$x[enr$cell_ID %in% a_ids]))
    expect_identical(p$labels$x, "mean neighborhood composition")
})

test_that("groups follow a custom order, or sort without clustering", {
    g <- .lp_gobject()
    p <- .lp(g, groups = "grp", labels = "lab",
        custom_group_order = c("c", "b", "a"))
    expect_identical(levels(p$data$group), c("c", "b", "a"))
    p2 <- .lp(g, groups = "grp", labels = "lab", cluster_groups = FALSE)
    expect_identical(levels(p2$data$group), c("a", "b", "c"))
})

test_that("a view counts only the cells it keeps", {
    g <- .lp_gobject()
    md <- GiottoClass::pDataDT(g)
    gv <- subset(g, lab == "x" | grp == "c", view = "v")
    kept <- md[lab == "x" | grp == "c"]$cell_ID
    got <- .lp(gv, groups = "grp", labels = "lab", view = "v")
    ref <- .lp(GiottoClass::subsetGiotto(g, cell_ids = kept),
        groups = "grp", labels = "lab")
    expect_equal(got$data[order(group, label), value],
        ref$data[order(group, label), value])
    expect_equal(.lp_value(got, "a", "x"), 1)
})

test_that("the dendrogram is drawn in the same plot, which takes layers", {
    skip_if_not_installed("ggdendro")
    g <- .lp_gobject()
    p <- plotLabelProportions(g, groups = "grp", labels = "lab",
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
    expect_s3_class(p, "ggplot")
    expect_true(any(vapply(p$layers, function(l)
        inherits(l$geom, "GeomSegment"), logical(1))))
    expect_identical(p$labels$fill, "lab")
    p2 <- p + ggplot2::labs(fill = "cell type")
    expect_identical(p2$labels$fill, "cell type")
    expect_no_error(ggplot2::ggplot_build(p2))
})

test_that("value_text_size sets the printed values", {
    g <- .lp_gobject()
    p <- .lp(g, groups = "grp", labels = "lab", value_text_size = 4)
    expect_identical(p$layers[[2]]$aes_params$size, 4)
})

test_that("inputs are checked", {
    g <- .lp_gobject()
    expect_error(.lp(g, groups = "grp"), "exactly one")
    expect_error(.lp(g, groups = "grp", labels = "lab",
        spat_enr_name = "nb"), "exactly one")
    expect_error(.lp(g, groups = "nope", labels = "lab"),
        "not in the cell metadata")
    expect_error(.lp(g, groups = "grp", labels = "lab",
        custom_group_order = c("a", "z")), "absent groups")
    tab <- GiottoClass::calculateLabelProportions(g, labels = "lab",
        group_method = "table", groups = "grp", output = "data.table",
        verbose = FALSE)
    expect_error(.lp(g, proportions = tab, value = "count"), "needs `labels`")
})
