test_that("scatterpie layer selects cell types by name, not position", {
    skip_if_not_installed("scatterpie")
    enr <- data.table::data.table(
        typeA = c(0.2, 0.5, 1),
        typeB = c(0.8, 0.5, 0),
        cell_ID = c("c1", "c2", "c3")
    )
    locs <- data.table::data.table(
        sdimx = c(1, 2, 3),
        sdimy = c(1, 2, 3),
        cell_ID = c("c1", "c2", "c3")
    )
    expect_no_error({
        pl <- plot_spat_scatterpie_layer_ggplot(
            ggobject = ggplot2::ggplot(),
            spatial_locations = locs,
            spatial_enrichment = enr
        )
        ggplot2::ggplot_build(pl)
    })
    fill_scale <- pl$scales$get_scales("fill")
    expect_setequal(names(fill_scale$palette(2)), c("typeA", "typeB"))
})
