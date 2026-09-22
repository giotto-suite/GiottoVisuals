# Tests for giottoMulti dispatch in GiottoVisuals (phase 4 of the
# gmulti-federation design — see GiottoClass/vignettes/DESIGN_gmulti_federation.md).
#
# Covers .resolve_samples (auto-injection + validation) and
# .gg_build_panel_child (joint cmeta projection via phase 3 access).

.mk_minimal <- function(ncell, nfeat) {
    m <- matrix(0, nrow = nfeat, ncol = ncell)
    rownames(m) <- paste0("f", seq_len(nfeat))
    colnames(m) <- paste0("c", seq_len(ncell))
    GiottoClass::createGiottoObject(expression = m, verbose = FALSE)
}

.mk_mg <- function() {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    GiottoClass::createGiottoMulti(list(B191 = g1, B215 = g2))
}


test_that(".resolve_samples defaults to all children", {
    mg <- .mk_mg()
    out <- GiottoVisuals:::.resolve_samples(mg, NULL, NULL, names(mg@objects))
    expect_setequal(out, c("B191", "B215"))
})

test_that(".resolve_samples respects an explicit subset", {
    mg <- .mk_mg()
    out <- GiottoVisuals:::.resolve_samples(mg, "B215", NULL,
        names(mg@objects))
    expect_identical(out, "B215")
})

test_that(".resolve_samples ':all:' sentinel expands to all children", {
    mg <- .mk_mg()
    out <- GiottoVisuals:::.resolve_samples(mg, ":all:", NULL,
        names(mg@objects))
    expect_setequal(out, c("B191", "B215"))
})

test_that(".resolve_samples errors clearly on unknown sample", {
    mg <- .mk_mg()
    expect_error(
        GiottoVisuals:::.resolve_samples(mg, "NOPE", NULL,
            names(mg@objects)),
        "not in giottoMulti"
    )
})

test_that(".gg_build_panel_child returns a giotto (not a giottoMulti)", {
    mg <- .mk_mg()
    panel <- GiottoVisuals:::.gg_build_panel_child(mg, "B191")
    expect_s4_class(panel, "giotto")
    expect_false(is(panel, "giottoMulti"))
})

test_that(".gg_build_panel_child has the child's own cells, not joint", {
    mg <- .mk_mg()
    panel <- GiottoVisuals:::.gg_build_panel_child(mg, "B191")
    panel_cm <- GiottoClass::getCellMetadata(panel, output = "data.table")
    expect_identical(nrow(panel_cm), 5L)  # B191 has 5 cells
    expect_true(all(panel_cm$cell_ID %in% paste0("c", 1:5)))
})

test_that(".gg_build_panel_child projects joint-only cmeta columns", {
    mg <- .mk_mg()
    # Inject a joint-only column ("leiden_clus") onto the gmulti's
    # joint cmeta — emulates what subsequent clustering would write.
    joint_cm <- GiottoClass::getCellMetadata(mg, output = "data.table")
    joint_cm[, leiden_clus := paste0("cl", rep(1:2, length.out = .N))]
    mg <- GiottoClass::setCellMetadata(mg,
        x = GiottoClass::createCellMetaObj(joint_cm))

    panel <- GiottoVisuals:::.gg_build_panel_child(mg, "B191")
    panel_cm <- GiottoClass::getCellMetadata(panel, output = "data.table")
    expect_true("leiden_clus" %in% names(panel_cm))
    expect_identical(nrow(panel_cm), 5L)
})

test_that(".gg_build_panel_child is a no-op when no joint-only columns exist", {
    mg <- .mk_mg()
    # No joint cmeta extras — panel cmeta = child's cmeta.
    panel <- GiottoVisuals:::.gg_build_panel_child(mg, "B191")
    panel_cols <- names(GiottoClass::getCellMetadata(panel, output = "data.table"))
    child_cols <- names(GiottoClass::getCellMetadata(mg@objects[["B191"]],
        output = "data.table"))
    expect_setequal(panel_cols, child_cols)
})

test_that(".gg_build_panel_child errors on unknown sample", {
    mg <- .mk_mg()
    expect_error(GiottoVisuals:::.gg_build_panel_child(mg, "NOPE"),
        "no child named")
})


# space -> panel set ###########################################################

# Which samples a `space` implies depends on whether its membership is
# closed. These lock that distinction down, because getting it wrong is
# silent both ways: reading membership off the wrong kind either drops
# panels a space covers, or draws an untransformed panel beside
# transformed ones.

.mk_mg_spaces <- function() {
    mg <- .mk_mg()
    # combined: B191 only, so its membership is a PROPER subset of the
    # children -- the case where deriving the panel set actually does
    # something.
    GiottoClass::giottoSpace(mg, "layout") <- GiottoClass::combinedSpace("B191")
    GiottoClass::giottoSpace(mg, "each") <- GiottoClass::spatShift(
        GiottoClass::perSampleSpace(), dx = 10)
    mg
}


test_that("a combinedSpace's membership is the default panel set", {
    mg <- .mk_mg_spaces()
    out <- GiottoVisuals:::.resolve_samples(mg, NULL, "layout",
        names(mg@objects))
    expect_identical(out, "B191")
})

test_that("':all:' under a combinedSpace still means the space's members", {
    mg <- .mk_mg_spaces()
    out <- GiottoVisuals:::.resolve_samples(mg, ":all:", "layout",
        names(mg@objects))
    expect_identical(out, "B191")
})

test_that("a sample outside a combinedSpace is refused", {
    mg <- .mk_mg_spaces()
    expect_error(
        GiottoVisuals:::.resolve_samples(mg, "B215", "layout",
            names(mg@objects)),
        "not in space"
    )
})

test_that("a perSampleSpace does not constrain the panel set", {
    # Its membership is open by design: an unscoped step applies to
    # whatever sample it meets, so `names()` on it lists what the recipe
    # mentions (here: nothing), not what it covers.
    mg <- .mk_mg_spaces()
    expect_length(names(GiottoClass::giottoSpace(mg, "each")), 0L)

    expect_setequal(
        GiottoVisuals:::.resolve_samples(mg, NULL, "each",
            names(mg@objects)),
        c("B191", "B215"))
    expect_identical(
        GiottoVisuals:::.resolve_samples(mg, "B215", "each",
            names(mg@objects)),
        "B215")
})

test_that("an unknown space name errors rather than plotting everything", {
    mg <- .mk_mg_spaces()
    expect_error(
        GiottoVisuals:::.resolve_samples(mg, NULL, "nope",
            names(mg@objects))
    )
})


# space -> panel contents ######################################################

# The panel set is only half of what `space =` means. These go through
# the real dispatcher because the other half -- that the frame reaches
# the panels at all -- is invisible to `.resolve_samples`: the recipe is
# applied once on the multi before the loop, since that is the last
# point where a child's sample name is still known.

.mk_mg_spatial <- function() {
    mk <- function(n, seed) {
        set.seed(seed)
        m <- matrix(stats::rpois(4L * n, 5), nrow = 4L)
        rownames(m) <- paste0("f", seq_len(4L))
        colnames(m) <- paste0("c", seq_len(n))
        GiottoClass::createGiottoObject(
            expression = m,
            spatial_locs = data.frame(sdimx = stats::runif(n),
                sdimy = stats::runif(n), cell_ID = colnames(m)),
            verbose = FALSE)
    }
    GiottoClass::createGiottoMulti(
        list(A = mk(20L, 1L), B = mk(15L, 2L), C = mk(12L, 3L)))
}


test_that("a combinedSpace narrows the panels actually drawn", {
    mg <- .mk_mg_spatial()
    GiottoClass::giottoSpace(mg, "layout") <- GiottoClass::spatShift(
        GiottoClass::combinedSpace(c("A", "B")), dx = 100, samples = "B")

    n_panels <- function(p) length(p$layers)  # one grob per sub-plot
    expect_identical(
        n_panels(spatPlot2D(mg, show_plot = FALSE, return_plot = TRUE,
            save_plot = FALSE)), 3L)
    expect_identical(
        n_panels(spatPlot2D(mg, space = "layout", show_plot = FALSE,
            return_plot = TRUE, save_plot = FALSE)), 2L)
})

test_that("a space reaches the panel, not just the panel set", {
    # A one-member space makes the composite a single ggplot, so the
    # coordinates it drew are readable off `$data`.
    mg <- .mk_mg_spatial()
    GiottoClass::giottoSpace(mg, "only_b") <- GiottoClass::spatShift(
        GiottoClass::combinedSpace("B"), dx = 100)

    native <- GiottoClass::getSpatialLocations(mg@objects$B,
        output = "data.table")
    p <- spatPlot2D(mg, space = "only_b", show_plot = FALSE,
        return_plot = TRUE, save_plot = FALSE)

    drawn <- p$layers[[1L]]$data
    expect_equal(sort(drawn$sdimx), sort(native$sdimx + 100))
})
