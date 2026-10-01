# The annotated cluster tree figure.
#
# It is three panels sharing one x axis, so the failures worth pinning are
# structural rather than cosmetic: it must render from a partially labelled
# answer, with or without a ladder, with or without evidence tracks.

.pct_gobject <- function(n_genes = 30L, n_cells = 60L, n_clus = 6L, seed = 7L) {
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
    list(gobject = g, n_clus = n_clus)
}

.pct_tree <- function(g) {
    Giotto::calculateClusterTree(g, cluster_column = "clus",
        expression_values = "raw", cor = "pearson", distance = "average")
}

.pct_labels <- function(tree) {
    ids <- as.character(utils::tail(seq_len(nrow(tree$merge)), 2L))
    list(
        clusters = stats::setNames(paste0("type_", tree$labels), tree$labels),
        nodes = stats::setNames(paste0("clade_", ids), ids)
    )
}


test_that("plotClusterTree renders with a ladder and tracks", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    tree <- .pct_tree(fx$gobject)
    gi <- Giotto::findGiniMarkers_one_vs_all(fx$gobject,
        cluster_column = "clus", expression_values = "raw",
        min_feats = 2, verbose = FALSE)

    p <- plotClusterTree(fx$gobject, cluster_column = "clus", tree = tree,
        labels = .pct_labels(tree), k = c(2L, 3L), gini_markers = gi,
        return_plot = TRUE, show_plot = FALSE, save_plot = FALSE)

    expect_s3_class(p, "ggplot")
})


test_that("plotClusterTree renders without a ladder or tracks", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    tree <- .pct_tree(fx$gobject)

    p <- plotClusterTree(fx$gobject, cluster_column = "clus", tree = tree,
        labels = .pct_labels(tree), return_plot = TRUE, show_plot = FALSE,
        save_plot = FALSE)
    expect_s3_class(p, "ggplot")
})


test_that("plotClusterTree renders with no node labels at all", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    tree <- .pct_tree(fx$gobject)
    labs <- list(clusters = stats::setNames(
        paste0("type_", tree$labels), tree$labels))

    p <- plotClusterTree(fx$gobject, cluster_column = "clus", tree = tree,
        labels = labs, k = c(2L, 4L), return_plot = TRUE, show_plot = FALSE,
        save_plot = FALSE)
    expect_s3_class(p, "ggplot")
})


test_that("plotClusterTree survives a sparsely labelled node set at every k", {
    # the ladder resolves a name per band the same way the annotation does, so
    # it has the same hazard: `[[` on a named vector errors for an absent name
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    tree <- .pct_tree(fx$gobject)
    labs <- list(
        clusters = stats::setNames(paste0("type_", tree$labels), tree$labels),
        nodes = stats::setNames("only_one", "1")
    )
    for (kk in 2:fx$n_clus) {
        expect_error(
            plotClusterTree(fx$gobject, cluster_column = "clus", tree = tree,
                labels = labs, k = kk, return_plot = TRUE, show_plot = FALSE,
                save_plot = FALSE),
            NA
        )
    }
})


test_that("no local variable shadows what [.data.table uses internally", {
    # `[.data.table` evaluates `setkeyv(ans, key)` when `by=` is given, and
    # resolves `key` against the calling frame. A local named `key` in
    # plotClusterTree made the grouped call fail trying to key its result by a
    # vector of cell type names. Pinned by name so it cannot come back.
    src <- deparse(GiottoVisuals::plotClusterTree)
    expect_false(any(grepl("^\\s*key\\s*<-", src)))
})


test_that("plotClusterTree refuses what it cannot draw", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    tree <- .pct_tree(fx$gobject)
    labs <- .pct_labels(tree)

    expect_error(
        plotClusterTree(fx$gobject, cluster_column = "clus", tree = list(1),
            labels = labs),
        "must be an `hclust`"
    )
    expect_error(
        plotClusterTree(fx$gobject, cluster_column = "clus", tree = tree,
            labels = labs, k = fx$n_clus + 5L),
        "must be between 1 and the number of clusters"
    )
})


test_that("plotClusterTree draws a giottoTree on its own", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    tree <- .pct_tree(fx$gobject)
    # no labels, no cluster_column: the tree records the column it was built on
    p <- plotClusterTree(fx$gobject, tree,
        return_plot = TRUE, show_plot = FALSE, save_plot = FALSE)
    expect_s3_class(p, "ggplot")
})


test_that("plotClusterTree needs cluster_column for a plain hclust", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    plain <- .pct_tree(fx$gobject)
    class(plain) <- "hclust"
    attr(plain, "params") <- NULL
    expect_error(plotClusterTree(fx$gobject, plain, show_plot = FALSE),
        "`cluster_column` is needed")
    p <- plotClusterTree(fx$gobject, plain, cluster_column = "clus",
        return_plot = TRUE, show_plot = FALSE, save_plot = FALSE)
    expect_s3_class(p, "ggplot")
})


test_that("plotClusterTree refuses a tree from another clustering", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    tree <- .pct_tree(fx$gobject)
    cx <- GiottoClass::getCellMetadata(fx$gobject, output = "cellMetaObj")
    cx[][clus == "6", "clus" := "7"]
    g2 <- GiottoClass::setCellMetadata(fx$gobject, cx, verbose = FALSE,
        initialize = FALSE)
    expect_error(plotClusterTree(g2, tree, show_plot = FALSE),
        "clusters with no leaf: 7.*leaves with no cells: 6")
})


test_that("plotClusterTree counts within the tree's view, and flags another", {
    skip_if_not_installed("ggdendro")
    fx <- .pct_gobject()
    g <- GiottoClass::addCellMetadata(fx$gobject, new_metadata = data.frame(
        cell_ID = GiottoClass::spatIDs(fx$gobject),
        keep = GiottoClass::pDataDT(fx$gobject)$clus %in% c("1", "2", "3")
    ))
    g <- subset(g, keep == TRUE, view = "v")
    g <- subset(g, clus != "none", view = "all")
    tree <- Giotto::calculateClusterTree(g, cluster_column = "clus",
        expression_values = "raw", distance = "average", view = "v")

    expect_message(
        p <- plotClusterTree(g, tree, return_plot = TRUE, show_plot = FALSE,
            save_plot = FALSE),
        "using view \"v\" recorded on the tree"
    )
    expect_s3_class(p, "ggplot")
    # the whole object brings back clusters the tree has no leaf for
    expect_error(expect_warning(
        plotClusterTree(g, tree, view = "all", show_plot = FALSE),
        "the tree was built with `view = \"v\"`"
    ), "clusters with no leaf")
})
