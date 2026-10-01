# showClusterHeatmap(tree = ): the tree supplies the branches and order on both
# axes; the correlations are computed here, from the cells and features the
# tree was built on unless the caller says otherwise.

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
    g <- GiottoClass::addCellMetadata(g,
        new_metadata = data.frame(cell_ID = colnames(m), clus = clus))
    list(gobject = g, mat = m, clus = clus)
}

.ch_feats <- paste0("g", c(1:12, 20:25))

.ch_tree <- function(g) {
    Giotto::calculateClusterTree(g, cluster_column = "clus",
        expression_values = "raw", feats = .ch_feats, distance = "average")
}

.ch_heat <- function(...) {
    showClusterHeatmap(..., show_plot = FALSE, return_plot = TRUE,
        save_plot = FALSE)
}


test_that("the tree's branches and order are drawn on both axes", {
    skip_if_not_installed("ComplexHeatmap")
    fx <- .ch_gobject()
    tree <- .ch_tree(fx$gobject)
    hm <- .ch_heat(fx$gobject, tree = tree)
    expect_identical(rownames(hm@matrix), tree$labels)
    expect_identical(colnames(hm@matrix), tree$labels)
    expect_identical(hm@row_dend_param$obj$merge, tree$merge)
    expect_identical(hm@column_dend_param$obj$merge, tree$merge)
})

test_that("colours are computed from the tree's features, not read off it", {
    skip_if_not_installed("ComplexHeatmap")
    fx <- .ch_gobject()
    tree <- .ch_tree(fx$gobject)
    hm <- .ch_heat(fx$gobject, tree = tree)

    # independent reference: per-cluster means over the tree's features
    sub <- as.matrix(fx$mat[.ch_feats, , drop = FALSE])
    means <- vapply(tree$labels, function(k) {
        rowMeans(sub[, fx$clus == k, drop = FALSE])
    }, numeric(nrow(sub)))
    expect_equal(hm@matrix, stats::cor(means), ignore_attr = TRUE)
})

test_that("other features override the tree's, with a warning", {
    skip_if_not_installed("ComplexHeatmap")
    fx <- .ch_gobject()
    tree <- .ch_tree(fx$gobject)
    expect_warning(
        hm <- .ch_heat(fx$gobject, tree = tree, feats = paste0("g", 1:30)),
        "the tree was built with `feats = "
    )
    # the same features in another order are not a conflict
    expect_no_warning(.ch_heat(fx$gobject, tree = tree,
        feats = rev(.ch_feats)))
})

test_that("a plain hclust needs cluster_column", {
    skip_if_not_installed("ComplexHeatmap")
    fx <- .ch_gobject()
    plain <- .ch_tree(fx$gobject)
    class(plain) <- "hclust"
    attr(plain, "params") <- NULL
    expect_error(.ch_heat(fx$gobject, tree = plain),
        "`cluster_column` is needed")
    hm <- .ch_heat(fx$gobject, tree = plain, cluster_column = "clus",
        expression_values = "raw")
    expect_identical(rownames(hm@matrix), plain$labels)
})

test_that("a tree from another clustering is refused", {
    skip_if_not_installed("ComplexHeatmap")
    fx <- .ch_gobject()
    tree <- .ch_tree(fx$gobject)
    cx <- GiottoClass::getCellMetadata(fx$gobject, output = "cellMetaObj")
    cx[][clus == "6", "clus" := "7"]
    g2 <- GiottoClass::setCellMetadata(fx$gobject, cx, verbose = FALSE,
        initialize = FALSE)
    expect_error(.ch_heat(g2, tree = tree),
        "clusters with no leaf: 7.*leaves with no cells: 6")
})

test_that("without a tree it still clusters the matrix itself", {
    skip_if_not_installed("ComplexHeatmap")
    fx <- .ch_gobject()
    hm <- .ch_heat(fx$gobject, cluster_column = "clus",
        expression_values = "raw")
    expect_s4_class(hm, "Heatmap")
    expect_false(inherits(hm@row_dend_param$obj, "giottoTree"))
})
