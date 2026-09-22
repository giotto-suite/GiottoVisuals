#' @include aux_visuals.R
NULL


# gmulti multi-sample dispatch for spatial plot functions
#
# Pattern: spatial plot functions that opt into multi-sample dispatch
# insert a guard at the top of their body, e.g.
#
#     spatPlot2D <- function(gobject, ..., view = NULL, space = NULL,
#                            samples = NULL) {
#         if (inherits(gobject, "giottoMulti")) {
#             return(.gg_multi_dispatch_spatial(
#                 plot_fn = spatPlot2D,
#                 named = mget(setdiff(names(formals()), "...")),
#                 dots = list(...),    # omit if the fn has no `...`
#                 gobject = gobject, view = view, space = space,
#                 samples = samples
#             ))
#         }
#         .gg_assert_giotto_single(gobject)
#         ...
#     }
#
# Three call-time axes (orthogonal, see DESIGN_gmulti_federation.md in
# GiottoClass for the full design):
#   - `view`   = subset/filter recipe (referenced by name from gobject@view)
#   - `space`  = coordinate-frame recipe (referenced by name from gobject@spaces);
#                applied to the panel children, and -- when its membership
#                is closed -- also the default panel set.
#   - `samples` = ad-hoc multi-sample selection (character vector, `":all:"`
#                sentinel, or NULL). NOT a slotted recipe; pure call-time arg.
#
# Resolution rules (samples auto-injection from space). Only a
# `combinedSpace` constrains the panel set; see `.resolve_samples` for why
# a `perSampleSpace` deliberately does not:
#   samples = NULL, space = NULL     -> all children, native frame
#   samples = NULL, combined "S"     -> S's members, in S's frame
#   samples = c(...), combined "S"   -> error if any sample is not in S
#   samples = NULL, perSample "S"    -> all children, each in S's frame
#   samples = c(...), space = NULL   -> those samples, native frame
#
# Per-panel call: each iteration extracts a single-sample "panel child"
# via `.gg_build_panel_child` — a giotto built from the gmulti's child
# augmented with joint slot projections (cmeta columns like leiden_clus
# from the joint level) sliced to that sample's contributions. Uses
# only exported GiottoClass APIs; no `:::` into internals.

#' @keywords internal
#' @noRd
.gg_multi_dispatch_spatial <- function(
        plot_fn, named, dots = list(), gobject, view, space, samples = NULL) {
    checkmate::assert_class(gobject, "giottoMulti")
    # A name, never an inline recipe -- see `.gg_materialize`.
    if (!is.null(view)) checkmate::assert_string(view, .var.name = "view")
    if (!is.null(space)) checkmate::assert_string(space, .var.name = "space")

    child_names <- names(gobject@objects)
    if (length(child_names) == 0L) {
        stop("[gmulti dispatch] giottoMulti has no child gobjects",
            call. = FALSE)
    }

    samples <- .resolve_samples(gobject, samples, space, child_names)

    # View and space both apply once across the gmulti before the panel
    # loop. The resolver's surviving-cell cache fills here, and each
    # per-child getter chain in the loop below pulls only the slice
    # belonging to its sample.
    #
    # `space` has to be applied HERE rather than forwarded to the panel
    # call: `materialize()` on a multi hands each child the recipe scoped
    # to its own name (`space_obj[samp]`), which is the only place the
    # sample identity needed to pick a per-sample step is still known. A
    # panel child is a plain `giotto` and has no name to resolve against.
    if (!is.null(view) || !is.null(space)) {
        gobject <- GiottoClass::materialize(gobject, view, space = space)
        # selectSamples in `view` may have narrowed children; reconcile
        # `samples` to the survivors. If a caller-supplied sample didn't
        # survive view narrowing, that's an error (silent drop is too
        # subtle).
        survivors <- names(gobject@objects)
        bad <- setdiff(samples, survivors)
        if (length(bad) > 0L) {
            stop(sprintf(paste(
                "[gmulti dispatch] sample(s) '%s' were filtered out by",
                "the applied view; cannot render panels for them.",
                "Surviving samples: %s",
                sep = " "
            ), paste(bad, collapse = ", "),
                paste(survivors, collapse = ", ")), call. = FALSE)
        }
    }

    # Title composition: user-supplied prefix + sample name as suffix.
    # Only applied when the plot function actually accepts a `title`
    # formal — some spatial plots (e.g. spatFeatPlot2D, spatCellPlot2D)
    # title sub-panels by feature name and don't expose a top-level
    # `title` arg.
    fn_formals <- names(formals(plot_fn))
    has_title <- "title" %in% fn_formals
    base_title <- if (has_title) named$title else NULL

    plots <- lapply(samples, function(s) {
        a <- named
        a$gobject <- .gg_build_panel_child(gobject, s)
        a$view <- NULL    # already applied
        a$space <- NULL   # already applied, per child, above
        a$samples <- NULL # consumed by dispatcher
        if (has_title) {
            a$title <- if (is.null(base_title) ||
                identical(base_title, "")) {
                s
            } else {
                paste(base_title, s)
            }
        }
        # Suppress per-panel show / save / return-override; the
        # dispatch returns one composite plot which the outer
        # plot_output_handler still processes.
        a$show_plot <- FALSE
        a$return_plot <- TRUE
        a$save_plot <- FALSE
        do.call(plot_fn, c(a, dots))
    })

    composite <- if (length(plots) == 1L) {
        plots[[1L]]
    } else {
        cowplot::plot_grid(
            plotlist = plots,
            ncol = set_default_cow_n_col(
                cow_n_col = named$cow_n_col,
                nr_plots = length(plots)
            ),
            rel_heights = named$cow_rel_h %null% 1,
            rel_widths = named$cow_rel_w %null% 1,
            align = named$cow_align %null% "h"
        )
    }

    # Route the composite through plot_output_handler so the user's
    # show_plot / save_plot / return_plot / save_param / default_save_name
    # are honored at the gmulti level. Per-panel values were already
    # suppressed in the lapply above.
    plot_output_handler(
        gobject = gobject,
        plot_object = composite,
        save_plot = named$save_plot,
        return_plot = named$return_plot,
        show_plot = named$show_plot,
        default_save_name = named$default_save_name %null%
            paste0(deparse(substitute(plot_fn)), "_gmulti"),
        save_param = named$save_param %null% list(),
        else_return = NULL
    )
}


# Resolve which samples the panel loop iterates over.
#
# - `samples = NULL` and no `space`        -> all children
# - `samples = NULL` and a combinedSpace   -> that space's members
# - `samples = NULL` and a perSampleSpace  -> all children
# - `samples = ":all:"`                    -> as the NULL cases above
# - `samples = c("A","B")`                 -> those samples; must exist on
#                                              the gmulti, and must be
#                                              members of a combinedSpace
#                                              if one was named
#
# All resolved sample names must be in `child_names` (i.e. real
# `@objects` keys). Defined-space names that overlap with sample names
# trigger no special handling here — the user is unambiguous because
# this arg is `samples`, not `space`.
#
# Why the two space kinds differ. A `combinedSpace` lays its members out
# relative to one another in one coordinate system, so its membership is
# CLOSED and declared: a sample outside the layout has no position in the
# frame, and drawing it would put an untransformed panel beside
# transformed ones. A `perSampleSpace` is OPEN by design — an unscoped
# step means "whatever sample it meets" — so `names()` on it lists the
# samples the recipe happens to mention, not the samples it covers.
# Deriving a panel set from that would silently drop children the space
# applies to perfectly well.
#
# This is narrower than §10.8's rule for artifact generators, where an
# explicit `samples =` may never escape the space at all. A panel set is
# not an artifact: nothing is persisted under a name that would later
# read as covering more than it does, and the scope is visible in the
# call that produced the figure. So the constraint is applied only where
# the frame itself makes an outside sample meaningless.
#' @keywords internal
#' @noRd
.resolve_samples <- function(gobject, samples, space, child_names) {
    # A named space that does not exist is an error, not a fallback to
    # "all children in the native frame" — `giottoSpace()` is the one
    # place a space name is checked and it diagnoses a name that is
    # really a sample or a group.
    space_samples <- NULL
    if (is.character(space) && length(space) == 1L) {
        space_obj <- GiottoClass::giottoSpace(gobject, space)
        if (inherits(space_obj, "combinedSpace")) {
            # `names()` on a space IS its membership. The `@samples` slot
            # this used to read no longer exists; the read sat inside a
            # `tryCatch` that swallowed the error, so the space silently
            # constrained nothing.
            space_samples <- names(space_obj)
            if (length(space_samples) == 0L) space_samples <- NULL
        }
    }

    if (is.null(samples)) {
        samples <- space_samples %||% ":all:"
    }
    if (identical(samples, ":all:")) {
        samples <- if (!is.null(space_samples)) space_samples else child_names
    }
    checkmate::assert_character(samples,
        min.len = 1L, any.missing = FALSE)

    # Sample names must be real children.
    bad <- setdiff(samples, child_names)
    if (length(bad) > 0L) {
        stop(sprintf(paste(
            "[gmulti dispatch] sample(s) '%s' not in giottoMulti.",
            "Available: %s",
            sep = " "
        ), paste(bad, collapse = ", "),
            paste(child_names, collapse = ", ")), call. = FALSE)
    }

    # If a defined space is in play, all samples must participate in it.
    if (!is.null(space_samples)) {
        unsupported <- setdiff(samples, space_samples)
        if (length(unsupported) > 0L) {
            stop(sprintf(paste(
                "[gmulti dispatch] sample(s) '%s' are not in space",
                "'%s' (its samples: %s). Drop them from `samples` or",
                "extend the space.",
                sep = " "
            ), paste(unsupported, collapse = ", "), space,
                paste(space_samples, collapse = ", ")), call. = FALSE)
        }
    }

    samples
}


# Build a single-giotto "panel child" for the per-panel plot fn call.
#
# Starts from `gobject@objects[[sample]]` (the actual child) and projects
# joint-level metadata columns onto it: when a column exists in the
# multi's joint @cell_metadata but not in the child's local cmeta, this
# adds it to the child via the public addCellMetadata API. Uses
# `getCellMetadata(mg, samples = ...)` from the access layer (phase 3) so
# only exported APIs are touched — no `:::` reach.
#
# The returned giotto is local to the panel iteration; mutations are
# not propagated back. Plot fns can call `getCellMetadata(g)` and see
# both child-local and joint columns uniformly.
#
# When the gmulti has no joint cmeta (or no joint-only columns), this
# returns the child unchanged.
#' @keywords internal
#' @noRd
.gg_build_panel_child <- function(gobject, sample) {
    child <- gobject@objects[[sample]]
    if (is.null(child)) {
        stop(sprintf("[gmulti dispatch] no child named '%s'", sample),
            call. = FALSE)
    }

    cell_ID <- NULL  # data.table NSE

    # Joint cmeta projection. Sliced to this sample only; cell_ID prefix
    # stripped so it lines up with child-local cell IDs.
    joint_cm <- tryCatch(
        GiottoClass::getCellMetadata(gobject,
            samples = sample, output = "data.table"),
        error = function(e) NULL)
    if (is.null(joint_cm) || nrow(joint_cm) == 0L) return(child)

    child_cm <- tryCatch(
        GiottoClass::getCellMetadata(child, output = "data.table"),
        error = function(e) NULL)
    if (is.null(child_cm)) return(child)

    joint_only <- setdiff(names(joint_cm), c(names(child_cm), "list_ID"))
    if (length(joint_only) == 0L) return(child)

    joint_cm <- data.table::copy(joint_cm)
    prefix <- paste0(sample, "::")
    joint_cm[, cell_ID := sub(paste0("^", prefix), "", cell_ID)]

    GiottoClass::addCellMetadata(child,
        new_metadata = joint_cm[, c("cell_ID", joint_only), with = FALSE],
        by_column = TRUE,
        column_cell_ID = "cell_ID")
}
