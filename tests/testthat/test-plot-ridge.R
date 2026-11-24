# ridgePlot(): one feature across groups draws a ridge per group

test_that("one feature with group_by gives a ridge per group in one plot", {
    skip_if_not_installed("GiottoData")
    skip_if_not_installed("ggridges")
    g <- GiottoData::loadGiottoMini("visium")
    pl <- ridgePlot(g, feats = "nr_feats", group_by = "leiden_clus",
        group_order = "mean", scale_axis = "identity",
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
    expect_s3_class(pl, "ggplot")
    md <- GiottoClass::pDataDT(g)
    means <- tapply(md$nr_feats, md$leiden_clus, mean)
    expect_identical(levels(pl$data$feat), names(sort(means)))
    expect_identical(pl$labels$y, "leiden_clus")
})

test_that("several features keep one panel per group", {
    skip_if_not_installed("GiottoData")
    skip_if_not_installed("ggridges")
    g <- GiottoData::loadGiottoMini("visium")
    pl <- ridgePlot(g, feats = c("Gna12", "Ccnd2"), group_by = "leiden_clus",
        group_by_subset = 1:2, show_plot = FALSE, return_plot = TRUE,
        save_plot = FALSE)
    expect_null(pl$data$feat) # a cowplot grid of per-group panels
    expect_error(ridgePlot(g, feats = c("Gna12", "Ccnd2"),
        group_by = "leiden_clus", group_as = "ridge", show_plot = FALSE,
        return_plot = TRUE, save_plot = FALSE), "one feature")
})
