#' @title showClusterDendrogram
#' @name showClusterDendrogram
#' @description Creates dendrogram for selected clusters.
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param expression_values expression values to use
#' (e.g. "normalized", "scaled", "custom")
#' @param cluster_column name of column to use for clusters
#' (e.g. "leiden_clus")
#' @param cor correlation score to calculate distance
#' (e.g. "pearson", "spearman")
#' @param distance distance method to use for hierarchical clustering,
#' default to "ward.D"
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param rotate rotate dendrogram 90 degrees
#' @inheritDotParams ggdendro::ggdendrogram
#' @details Expression correlation dendrogram for selected clusters.
#' @returns ggplot
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#' showClusterDendrogram(g, cluster_column = "leiden_clus")
#'
#' @export
showClusterDendrogram <- function(gobject,
    spat_unit = NULL,
    feat_type = NULL,
    expression_values = c("normalized", "scaled", "custom"),
    cluster_column,
    cor = c("pearson", "spearman"),
    distance = "ward.D",
    h = NULL,
    h_color = "red",
    rotate = FALSE,
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    default_save_name = "showClusterDendrogram",
    view = NULL,
    space = NULL,
    ...) {
    # verify if optional package is installed
    package_check(pkg_name = "ggdendro", repository = "CRAN")
    gobject <- .gg_materialize(gobject, view, space,
        slots = c("cell_metadata", "expression", "spatial_enrichment"))

    values <- match.arg(
        expression_values,
        unique(c(
            "normalized", "scaled", "custom",
            expression_values
        ))
    )

    # Set feat_type and spat_unit
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type
    )

    metatable <- calculateMetaTable(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        expression_values = values,
        metadata_cols = cluster_column
    )

    pl <- create_cluster_dendrogram(
        data = metatable,
        clus_col = "uniq_ID",
        var_col = "variable",
        val_col = "value",
        cor = cor,
        distance = distance,
        h = h,
        h_color = h_color,
        rotate = rotate,
        ...
    )

    return(plot_output_handler(
        gobject = gobject,
        plot_object = pl,
        save_plot = save_plot,
        return_plot = return_plot,
        show_plot = show_plot,
        default_save_name = default_save_name,
        save_param = save_param,
        else_return = NULL
    ))
}


#' @name create_cluster_dendrogram
#' @title Create clustered expression dendrogram
#' @description Create a dendrogram based on a data.table with columns for
#' cluster ID, variables, and their values. If no specific values are provided
#' for the 'col' params then they will be assumed as 1. clus_col, 2. var_col,
#' 3. val_col
#' @param data data.table. Should include columns with clustering labels,
#' variables that are being clustered, and the values of those clusters.
#' @param clus_col character. name of column with clustering info
#' @param var_col character. name of column with variable name
#' @param val_col character. name of column with values info
#' @param cor correlation score to calculate distance
#' (e.g. "pearson", "spearman")
#' @param distance distance method to use for hierarchical clustering,
#' default to "ward.D"
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param rotate rotate dendrogram 90 degrees
#' @inheritDotParams ggdendro::ggdendrogram
#' @returns ggdendrogram
#' @examples
#' g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
#'
#' g_expression <- head(GiottoClass::getExpression(g, output = "matrix"))
#' g_expression_df <- as.data.frame(as.matrix(g_expression))
#' g_expression_df$feat_ID <- rownames(g_expression)
#'
#' g_expression_melt <- data.table::melt(g_expression_df,
#'     id.vars = "feat_ID",
#'     measure.vars = colnames(g_expression), variable.name = "cell_ID",
#'     value.name = "raw_expression"
#' )
#'
#' create_cluster_dendrogram(data.table::as.data.table(g_expression_melt),
#'     var_col = "cell_ID", clus_col = "feat_ID", "raw_expression"
#' )
#'
#' @export
create_cluster_dendrogram <- function(
        data,
        clus_col = names(data)[[1]],
        var_col = names(data)[[2]],
        val_col = names(data)[[3]],
        cor = c("pearson", "spearman"),
        distance = "ward.D",
        h = NULL,
        h_color = "red",
        rotate = FALSE,
        ...) {
    checkmate::assert_data_table(data)
    checkmate::assert_character(clus_col)
    checkmate::assert_character(var_col)
    checkmate::assert_character(val_col)
    cor <- match.arg(cor, c("pearson", "spearman"))

    dcast_metatable <- data.table::dcast.data.table(
        data = data,
        formula = paste0(var_col, "~", clus_col),
        value.var = val_col
    )
    testmatrix <- dt_to_matrix(x = dcast_metatable)

    # correlation
    cormatrix <- cor_flex(x = testmatrix, method = cor)
    cordist <- stats::as.dist(1 - cormatrix, diag = TRUE, upper = TRUE)
    corclus <- stats::hclust(d = cordist, method = distance)

    cordend <- stats::as.dendrogram(object = corclus)

    # plot dendrogram
    pl <- ggdendro::ggdendrogram(data = cordend, rotate = rotate, ...)

    # add horizontal or vertical lines
    if (!is.null(h)) {
        pl <- pl + ggplot2::geom_hline(yintercept = h, col = h_color)
    }

    pl
}
