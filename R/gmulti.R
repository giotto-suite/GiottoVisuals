#' @include aux_visuals.R
NULL


# gmulti multi-sample dispatch for spatial plot functions
#
# Pattern: spatial plot functions that opt into multi-sample dispatch
# insert a guard at the top of their body, e.g.
#
#     spatPlot2D <- function(gobject, ..., view = NULL, space = NULL) {
#         if (inherits(gobject, "giottoMulti")) {
#             return(.gg_multi_dispatch_spatial(
#                 plot_fn = spatPlot2D,
#                 named = mget(setdiff(names(formals()), "...")),
#                 dots = list(...),    # omit if the fn has no `...`
#                 gobject = gobject, view = view, space = space
#             ))
#         }
#         .gg_assert_giotto_single(gobject)
#         ...
#     }
#
# The helper extracts each requested child from `gobject@objects` and
# calls `plot_fn` on it as a single `giotto` — getters work without any
# space knowledge in the plot function itself. Panels are combined via
# `cowplot::plot_grid` reusing the plot function's existing
# `cow_n_col` / `cow_rel_h` / `cow_rel_w` / `cow_align` args.
#
# v1 scope: per-child dispatch only. Defined-space cross-sample panels
# (rendering N children together in one panel under a shared coord
# frame) are deferred until gmulti-side getters can flatten subobject
# lists into single subobjects. Until then, those entries error with
# a clear message.

#' @keywords internal
#' @noRd
.gg_multi_dispatch_spatial <- function(
        plot_fn, named, dots = list(), gobject, view, space) {
    checkmate::assert_class(gobject, "giottoMulti")

    child_names <- names(gobject@objects)
    if (length(child_names) == 0L) {
        stop("[gmulti dispatch] giottoMulti has no child gobjects",
            call. = FALSE)
    }

    # Default + `:all:` expansion → vector of all child names.
    if (is.null(space)) space <- ":all:"
    if (identical(space, ":all:")) space <- child_names
    checkmate::assert_character(space,
        min.len = 1L, any.missing = FALSE)

    # v1: only sample-name entries are supported. Surface a clear
    # error if the user passes a defined-space name (which is a real
    # thing on the gmulti, just not yet a renderable panel target).
    space_names <- tryCatch(
        GiottoClass::giottoSpaces(gobject),
        error = function(e) character()
    )
    bad <- setdiff(space, child_names)
    bad_space <- intersect(bad, space_names)
    if (length(bad_space) > 0L) {
        stop(sprintf(paste(
            "[gmulti dispatch] '%s' is a defined giottoSpace.",
            "Cross-sample defined-space panels are not yet supported",
            "(v1 dispatches per child only).",
            "Use sample names or `:all:` for now.",
            sep = " "
        ), bad_space[[1L]]), call. = FALSE)
    }
    if (length(bad) > 0L) {
        stop(sprintf(paste(
            "[gmulti dispatch] '%s' is not a sample name in this",
            "giottoMulti. Available samples: %s",
            sep = " "
        ), bad[[1L]], paste(child_names, collapse = ", ")),
        call. = FALSE)
    }

    # View applies once across the gmulti before the panel loop. The
    # resolver's surviving-cell cache fills here, and each per-child
    # getter chain in the loop below pulls only the slice belonging
    # to its sample.
    if (!is.null(view)) {
        gobject <- GiottoClass::materialize(gobject, view, space = NULL)
    }

    # Title composition: user-supplied prefix + sample name as suffix.
    # Only applied when the plot function actually accepts a `title`
    # formal — some spatial plots (e.g. spatFeatPlot2D, spatCellPlot2D)
    # title sub-panels by feature name and don't expose a top-level
    # `title` arg.
    fn_formals <- names(formals(plot_fn))
    has_title <- "title" %in% fn_formals
    base_title <- if (has_title) named$title else NULL

    plots <- lapply(space, function(child) {
        a <- named
        a$gobject <- gobject@objects[[child]]
        a$view <- NULL    # already applied
        a$space <- NULL   # single-sample below
        if (has_title) {
            a$title <- if (is.null(base_title) ||
                identical(base_title, "")) {
                child
            } else {
                paste(base_title, child)
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

    if (length(plots) == 1L) return(plots[[1L]])

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
