
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GiottoVisuals <img src="man/figures/logo.png" align="right" alt="" width="160" />

<!-- badges: start -->

![Version](https://img.shields.io/github/r-package/v/giotto-suite/GiottoVisuals)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![codecov](https://codecov.io/github/giotto-suite/GiottoVisuals/graph/badge.svg?token=F4ga1Ahbhw)](https://codecov.io/github/giotto-suite/GiottoVisuals)
[![R-CMD-check](https://github.com/giotto-suite/GiottoVisuals/actions/workflows/main_check.yml/badge.svg)](https://github.com/giotto-suite/GiottoVisuals/actions/workflows/main_check.yml)
[![GitHub
issues](https://img.shields.io/github/issues/giotto-suite/Giotto)](https://github.com/giotto-suite/Giotto/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/giotto-suite/GiottoVisuals)](https://github.com/giotto-suite/GiottoVisuals/pulls)
<!-- badges: end -->

GiottoVisuals contains the main plotting functions of Giotto Suite

## Installation

You can install GiottoVisuals like:

``` r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak("giotto-suite/GiottoVisuals")
```

## Script Organization by Prefixes:

`aux_` - auxilliary and meta functionality  
- `aux_output.R` - plot output handling  
- `aux_save.R` - plot saving  
- `aux_defaults.R` - set plotting defaults  
- `aux_visuals.R` - general auxilliary plotting functions

`plot_` - general plotting functions organized by type of plot  
- `plot_scatter.R` - scatter plots  
- `plot_violin.R` - violin plots  
- `plot_heatmap.R` - heatmaps  
- `plot_dendrogram.R` - dendrograms

`vis_` - specific plotting functions organized by topic - `vis_hvf.R` -
highly variable features plots  
- `vis_pc.R` - principle components plots  
- `vis_spatial.R` - - `vis_spatial_in_situ.R` -

`gg_` - ggplot2 related  
- `gg_info_layers.R` - modular internals for plotting the different
layers of information  
- \`\`

Other:  
- `color_palettes.R` - color palettes and factory functions


Package:  
- `package_imports.R` - all imports  
- `globals.R` - global variables  
- `dd.R` - dummy documentation for inheriting  
- `zzz.R` - onloads
