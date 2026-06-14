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
