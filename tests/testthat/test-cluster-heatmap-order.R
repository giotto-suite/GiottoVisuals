# showClusterHeatmap(cluster_custom_order = ): the clusters are laid out in the
# given order and not clustered, so no dendrogram is drawn -- the order is an
# arrangement, and the plot makes no claim about structure.

.ch_gobject <- function(n_genes = 30L, n_cells = 90L, n_clus = 6L, seed = 5L) {
    set.seed(seed)
    m <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.7,
        rand.x = function(n) as.double(rpois(n, 4L) + 1L))
    rownames(m) <- paste0("g", seq_len(n_genes))
    colnames(m) <- paste0("c", seq_len(n_cells))
    clus <- as.character(rep(seq_len(n_clus), length.out = n_cells))
    for (k in seq_len(n_clus)) {
        gi <- ((k - 1L) * 3L + 1L):(k * 3L)
        m[gi, clus == as.character(k)] <- m[gi, clus == as.character(k)] + 20
    }
    g <- GiottoClass::createGiottoObject(expression = m)
    GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(cell_ID = colnames(m), clus = clus))
}

.ch_heat <- function(...) {
    showClusterHeatmap(..., cluster_column = "clus",
        expression_values = "raw",
        show_plot = FALSE, return_plot = TRUE, save_plot = FALSE)
}


test_that("the clusters follow the given order, with no dendrogram", {
    skip_if_not_installed("ComplexHeatmap")
    g <- .ch_gobject()
    ord <- c("4", "1", "6", "2", "5", "3")
    hm <- .ch_heat(g, cluster_custom_order = ord)
    expect_identical(rownames(hm@matrix), ord)
    expect_identical(colnames(hm@matrix), ord)
    expect_false(hm@row_dend_param$cluster)
    expect_false(hm@column_dend_param$cluster)
})

test_that("the order only rearranges the correlations", {
    skip_if_not_installed("ComplexHeatmap")
    g <- .ch_gobject()
    ord <- c("4", "1", "6", "2", "5", "3")
    plain <- .ch_heat(g)
    ordered <- .ch_heat(g, cluster_custom_order = ord)
    expect_equal(ordered@matrix, plain@matrix[ord, ord])
    # and without an order it still clusters, drawing its own dendrogram
    expect_s3_class(plain@row_dend_param$obj, "hclust")
})

test_that("an order must name every cluster shown, once", {
    skip_if_not_installed("ComplexHeatmap")
    g <- .ch_gobject()
    expect_error(.ch_heat(g, cluster_custom_order = c("1", "2", "3")),
        "missing: 4, 5, 6")
    expect_error(.ch_heat(g, cluster_custom_order = c(1:6, 9)),
        "not in the data: 9")
    expect_error(.ch_heat(g, cluster_custom_order = c(1:6, 1)),
        "every cluster shown, once")
})
