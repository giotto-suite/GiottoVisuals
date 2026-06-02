# Giotto plot output handler #
# -------------------------- #


#' @name plot_output_handler
#' @title Plotting output handler
#' @description
#' Simple wrapper for handling Giotto's framework for plotting outputs and
#' saving Plotting functions should return using this handler.
#' @param gobject giotto object
#' @param plot_object plot object
#' @param save_plot logical. (defaults to instructions setting) whether to save
#' plot
#' @param return_plot logical. (defaults to instructions setting) whether to
#' return the plot object
#' @param show_plot logical. (defaults to instructions setting) whether to print
#' the plot in display device
#' @param default_save_name default name to use (should be set by the specific
#' plotting function. Users should not interact with this parameter)
#' @param save_param list of plot saving parameters see [showSaveParameters()]
#' @param else_return optional. What should be returned instead if the plot
#' object is not returned
#' @keywords internal
#' @returns plot object
#' @examples
#' g <- GiottoData::loadGiottoMini("vizgen")
#' g_spatplot <- spatPlot2D(g, return_plot = TRUE)
#'
#' plot_output_handler(g, plot_object = g_spatplot, save_plot = FALSE)
#'
#' @export
plot_output_handler <- function(gobject,
    plot_object,
    save_plot = NULL,
    return_plot = NULL,
    show_plot = NULL,
    default_save_name = NULL,
    save_param = list(),
    else_return = NULL) {
    # Accepts giotto or giottoMulti — instructions() works on both
    # via the gAny base. Plot fns that route their data through
    # combine* / spatValues / gmulti-aware getters can hand back a
    # gmulti here without per-sample splitting.
    checkmate::assert_class(gobject, "gAny")

    ## output settings detection ##
    show_plot <- show_plot %null% instructions(gobject, param = "show_plot")
    save_plot <- save_plot %null% instructions(gobject, param = "save_plot")
    return_plot <- return_plot %null%
        instructions(gobject, param = "return_plot")

    ## handle outputs --------------------------------------------------- ##

    ## print plot ##
    if (show_plot) {
        print(plot_object)
    }

    ## save plot ##
    if (save_plot) {
        checkmate::assert_character(default_save_name)
        checkmate::assert_list(save_param)

        data_param <- list(
            gobject = gobject,
            plot_object = plot_object,
            default_save_name = default_save_name
        )

        do.call("all_plots_save_function", args = c(data_param, save_param))
    }

    ## return plot ##
    if (return_plot) {
        return(invisible(plot_object))
    } else {
        return(invisible(else_return))
    }
}


# TODO split into two functions
# plot_output_handler_read = function(
#         gobject,
#         plot_object,
#         save_plot = NA,
#         return_plot = NA,
#         show_plot = NA,
#         default_save_name,
#         save_param,
#         else_return = NULL) {
#     checkmate::assert_class(gobject, 'giotto')
#     checkmate::assert_character(default_save_name)
#     checkmate::assert_list(save_param)
#
#     instr = instructions(gobject)
#
#     ## output settings detection ##
#     # IF setting is NA then the appropriate setting from gobject instructions
#     # will be checked and used.
#     # IF setting is NOT NA then the provided value will be used directly.
#     show_plot = ifelse(is.na(show_plot),
#                        readGiottoInstructions(instr, param = 'show_plot'),
#.                       show_plot)
#     save_plot = ifelse(is.na(save_plot),
#                        readGiottoInstructions(instr, param = 'save_plot'),
#.                       save_plot)
#     return_plot = ifelse(is.na(return_plot),
#                          readGiottoInstructions(instr, param = 'return_plot'),
#                          return_plot)
# }


# plot_output_handler_do = function(gplot_out) {
#     checkmate::assert_class(gplot_out, 'giottoPlotOutput')
#
#     ## print plot ##
#     if(gplot_out$show_plot) {
#         print(gplot_out$plot_object)
#     }
#
#     ## save plot ##
#     if(gplot_out$save_plot) {
#         do.call('all_plots_save_function',
#                 c(list(gobject = instr,
#                        plot_object = gplot_out$plot_object,
#                        default_save_name = default_save_name),
#                   save_param)
#         )
#     }
#
#     ## return plot ##
#     if(return_plot) {
#         invisible(return(plot_object))
#     } else {
#         return(else_return)
#     }
# }
