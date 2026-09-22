# GiottoVisuals 0.2.16 (in development)

## enhancements
- plot functions accept `view` and `space`: the name of a view or coordinate-frame recipe slotted on the object. The object is pre-narrowed through `GiottoClass::materialize()` once, so every slot read in a plot shares one resolver pass. Added to `spatPlot2D()`, `spatInSituPlotPoints()`, `spatInSituPlotHex()`, `spatInSituPlotDensity()`, `spatDeconvPlot()`, `spatFeatPlot2D()`, `spatCellPlot2D()` and the heatmap / violin / ridge / dot / dendrogram family.
- spatial plot functions take a `giottoMulti` and draw one panel per sample, composed with `cowplot::plot_grid()`. `samples` selects which children to draw; `":all:"` is the explicit form of the default.
- non-spatial and dim-reduction plots (`dotPlot()`, `violinPlot()`, `ridgePlot()`, `plotHeatmap()`, `plotMetaDataHeatmap()`, `showClusterHeatmap()`, `showClusterDendrogram()`, `dimPlot2D()`, `plotUMAP()`, `plotTSNE()`, `plotPCA()`) accept a `giottoMulti` directly, pooling cells across samples. Their getters already returned joint subobjects; only the class guard was rejecting them.
- `plot_output_handler()` accepts a `giottoMulti` as well as a `giotto`.

## bug fixes
- `show_network = TRUE` draws again. A spatial network has been stored as an igraph since GiottoClass 0.6.0, so `getSpatialNetwork(output = "networkDT")` returns an edge list with no coordinates, while the drawing code still expected `sdimx_begin` / `sdimy_begin` / `sdimx_end` / `sdimy_end`. The endpoints are now attached by `GiottoClass::annotateSpatialNetwork()` at draw time. Affected `spatPlot2D()`, `spatFeatPlot2D()` and the plotly 2D/3D family. The failure was at render rather than at fetch, so a plot object built without error and only broke when printed or saved.

## changes
- minimum GiottoClass version is `>= 0.7.2`, for `annotateSpatialNetwork()`'s optional annotations.
- a network drawn under a `view` shows only edges whose endpoints are drawn. This falls out of the endpoint lookup above — the view narrows the spatial locations and the edges follow — and matches what `select_cells` has always done. A `space` needs no handling: the locations the edges join to are already transformed.
- a `space` named on a spatial plot over a `giottoMulti` sets the panel set only when its membership is closed. A `combinedSpace` lays its members out relative to one another, so it defines the panels and refuses a `samples` outside it; a `perSampleSpace` is open by design, so it transforms every panel without constraining which are drawn.

# GiottoVisuals 0.2.15 (2026/05/14)

## bug fixes
- fix color gradient error when GiottoVisuals is not attached (loaded as dependency or via `::`): `giotto.color_cd_pal` and `giotto.color_cs_pal` options now have inline fallback defaults in `set_default_color_continuous()`
- `aes_string2()` fixed handling for column names that parse as numeric literals (e.g. "1", "2") that caused plotting to fail

## changes
- `auto_image_resample()` rewritten to use `terra::window()` instead of the two-method crop/oversample approach, simplifying the implementation and avoiding materialization of large crops to disk. `giottoAffineImage` is now also supported. Removed params: `flex_resample`, `max_crop`, `max_resample_scale` and their corresponding global options `giotto.plot_img_max_crop` and `giotto.plot_img_max_resample_scale`.
- feature value collection in `spatFeatPlot2D_single()`, `dimFeatPlot2D()`, `spatFeatPlot3D()`, `dimFeatPlot3D()`, `spatDimFeatPlot3D()`, and `violinPlot()` now uses `spatValues()`, replacing manual expression matrix extraction and transposition
- minimum GiottoClass version bumped to `>= 0.5.1`

# GiottoVisuals 0.2.14 (2025/11/19)

## bug fixes
- allow `coord_fix_ratio = TRUE` usage again after ggplot2 changes

## changes
- change `spatInSituPlotPoints()` param `spat_enr_names` to `spat_enr_name`

## enhancements
- `show_axes` param for `spatInSituPlotPoints()`

# GiottoVisuals 0.2.13 (2025/09/30)

## bug fixes
- update for removal of `ggplot2::aes_string()` in ggplot2 v4 [#1192](https://github.com/drieslab/Giotto/issues/1192) by jdub0211

## changes
- replacement of argument radius by pie_scale in `spatDeconvPlot()` to modify pie size.

# GiottoVisuals 0.2.12 (2025/05/06)

## new
- `ridgePlot()`

## bug fixes
- fix `spatInSituPlotDensity()` `polygon_alpha` default setting when polygons plotted last

## enhancements
- requested images not within the plotting window are now handled as warning instead of erroring
- `y_axis_title` param added to `plotMetaDataCellsHeatmap()`
- `spatInSituPlotDensity()` `xlim`/`ylim` now crops on polys
- `spatInSituPlotDensity()` feature plotting via `spatValues()`retrieval

# GiottoVisuals 0.2.11 (2025/01/17)

## enhancements
- improve efficiency of `spatInSituPlotDensity()`
- `use_overlap` param for `spatInSituPlotDensity()` to allow the unoverlapped data to be used for density

## bug fixes
- fix color code/name checking to allow hex codes in `cell_color` plotting param.

# GiottoVisuals 0.2.10 (2024/12/04)

## bug fixes
- fix missing `e` extent error during `spatInSituPlotPoints()` with `show_image = TRUE` introduced in v0.2.9 [#96](https://github.com/drieslab/GiottoVisuals/issues/96) by kunleng

# GiottoVisuals 0.2.9 (2024/11/18)

## bug fixes
- `spatInSituPlotPoints()` `xlim` and `ylim` params now also affect image resampling

## enhancements
- `group_by` and `group_by_subset` for `dimFeatPlot2D()` [#1069](https://github.com/drieslab/Giotto/issues/1069) by xhNorthwestern
- image plotting will now scale values to the 99th percentile when a specific max intensity is not set.

# GiottoVisuals 0.2.8 (2024/11/14)

## enhancements
- `spatPlot3D()` now works with non-categorical data. Uses same params as `spatPlot2D()` for color/gradient selection

## changes
- `spatPlot3D()` `other_point_size` default changed to 3 from 0.5 which was too hard to notice.


# GiottoVisuals 0.2.7 (2024/11/08)

## enhancements
- scattermore plotting for `spatPlot2D()` and `spatFeatPlot2D()`

## bug fixes
- fix color scaling bug for `spatFeatPlot2D()` introduced in 0.2.6

# GiottoVisuals 0.2.6 (2024/10/27)

## new
- `combine_aes()` for combining ggplot2 mappings in a last item wins manner. Allows more than 2 aes to be combined at the same time and provides an optional informative warning
- `gg_params()` for combining and sorting ggplot2 mappings, toplevel params, and data into an arg list for plotting.


# GiottoVisuals 0.2.5 (2024/09/27)

## new
- `gpsparam()` for generating an object containing plot saving parameters
- `cowplot::plot_grid()` re-exported
- `dotPlot()` visualization

## bug fixes
- fix `coord_fix_ratio` passing in `spatFeatPlot2D()` [#972](https://github.com/drieslab/Giotto/issues/972) by kaizen89

# GiottoVisuals 0.2.4 (2024/07/26)

## enhancements
- `giottoLargeImage` `max_window` and `colors` slot info is now followed during ggplot plotting
- `giottoAffineImage` compatibility for giotto ggplot2 plotting functions
- `gg_annotation_raster()` now also performs `geom_blank()` with the extent provided through `ext` param. This can be turned off by setting `geom_blank() = FALSE`

## new
- `geom_text_repel()` and `geom_label_repel()` from `ggplot2` are now re-exported

# GiottoVisuals 0.2.3 (2024/05/28)

## bug fixes
- hotfix for sprintf numerics usage %d -> %f

# GiottoVisuals 0.2.2 (2024/05/22)

## breaking changes
- option `"giotto.color_c_pal"` has been split into `"giotto.color_cd_pal"` and `"giotto.color_cs_pal"` for divergent and sequential continuous palettes

## enhancements
- some options are now initialized on library load to make them easier to access

# GiottoVisuals 0.2.0 (2024/05/13)

## changes
- _GiottoClass >= 0.3.0_ is now required. (object no longer has `@largeImages` slot)
- `largeImage_name` arg is deprecated. `image_name` arg should be used instead
- spatial y dim values are not vertically plotted
- deprecated `spatDimGenePlot3D()`, `spatGenePlot3D()`, `dimGenePlot3D` in favor of the 'feat' naming scheme

## bug fixes
- removed centroid point that unintentionally gets plotted in polygons

## enhancements
- `theme_param` arg has been added to some plotting functions that provides access to `ggplot2::theme()` finetuning
- `gg_annotation_raster()` now performs plot extent detection and resampling of largeImages
- resampling args are now globally settable (see `?auto_image_resample`)
- `polygon_alpha` default is now 1 when plotting polys without images, and 0.5 when plotting with images.


# GiottoVisuals 0.1.7 (2024/03/11)

## Changes
- Change default values show_plot = NA, return_plot = NA, save_plot = NA -> show_plot = NULL, return_plot = NULL, save_plot = NULL.
- Add examples to exported functions
- Update License

# GiottoVisuals 0.1.6 (2024/03/07)

## bug fixes
- fix `spatPlot2D()` and `spatFeatPlot2D()` largeImage plotting when using `group_by`
- fix `edge_alpha` arg in `spatFeatPlot2D()`

## new
- `mixRGB()` vectorized additive mixing in RGB space
- `mixHSV()` vectorized color mixing in HSV space



# GiottoVisuals 0.1.4 (2024/01/25)
## bug fixes
- fix plotting color gradient when using arg `point_shape = "no_border"`
- fix image NA value [#865](https://github.com/drieslab/Giotto/issues/865) by rbutleriii

## new
- `gg_annotation_raster()` internal generic for adding a `giottoImage`, `giottoLargeImage`, or `list` thereof to a `ggplot` object

# GiottoVisuals 0.1.2 (2023/01/12)

## bug fixes
- fix `ComplexHeatmap` saving size issue [#861](https://github.com/drieslab/Giotto/issues/861) by rbutleriii

# GiottoVisuals 0.1.1 (2023/12/16)

## Breaking Changes
- Removed: `getDistinctColors()` to *GiottoUtils*

## Added
- Add: `getDistinctColors()` and `getRainbowColors()` as re-exports from *GiottoUtils*

## Changes
- Changed: Package internal functions now have `.` prefix
- Changed: *GiottoUtils* and *GiottoClass* moved to Imports
- Changed: Other internal function naming changes to bring in line with updates to *GiottoUtils*


# GiottoVisuals 0.1.0 (23/11/29)

Initial release
