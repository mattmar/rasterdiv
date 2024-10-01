
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![cran
version](http://www.r-pkg.org/badges/version/rasterdiv)](https://cran.r-project.org/package=rasterdiv)
[![Downloads
(monthly)](https://cranlogs.r-pkg.org/badges/last-month/rasterdiv?color=brightgreen)](https://cranlogs.r-pkg.org/badges/last-month/rasterdiv)
[![Downloads
(total)](https://cranlogs.r-pkg.org/badges/grand-total/rasterdiv?color=brightgreen)](https://www.r-pkg.org/pkg/rasterdiv)
<!-- badges: end -->

# rasterdiv

## Installation

### Stable version from [CRAN](https://cran.microsoft.com/):

``` r
install.packages("rasterdiv")
```

### Development version from [GitHub](https://github.com/):

``` r
# install.packages("remotes")
remotes::install_github("mattmar/rasterdiv")
```

## What is rasterdiv?

***rasterdiv*** is an R package that provides functions to apply
diversity indexes based on Information Theory on RasterLayer or
numerical matrices. Supported indexes include:

  - Parametric Rao’s quadratic entropy (classical and multidimensional);
  - Shannon’s diversity index;
  - Pielou’s evenness index;
  - Hill’s generalised entropy;
  - Rényi’s generalised entropy;
  - Berger-Parker index;
  - Cumulative Residual Entropy (CRE)

## New Features

  - Time Series Integration: The package recently integrated time series
    of remote sensing data to derive Rao’s index, for example
    considering plant phenology. See the vignette: [Multidimension Rao’s
    Index](https://mattmar.github.io/rasterdiv/articles/rasterdiv_03_Advanced_multidimension_Rao.html).
  - Helical Graphs: Prepare and plot time series with “helical graphs”.
    See the vignette: [Visualising rasterdiv indexes with Helical
    Plots](https://mattmar.github.io/rasterdiv/articles/rasterdiv_05_Helical_Plots.html).

## Key Features

  - Moving Window and Area-Based Calculations: Derive indexes using a
    “moving window” approach (all functions) or area-based Parametric
    Rao’s quadratic entropy via
    [*paRao()*](https://mattmar.github.io/rasterdiv/articles/rasterdiv_area_based_Rao.html)

## Learn More on **rasterdiv**

  - [***rasterdiv**—An Information Theory tailored R package for
    measuring ecosystem heterogeneity from space: To the origin and
    back*](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13583);
  - [*Measuring diversity from space: a global view of the free and open
    source rasterdiv R package under a coding
    perspective*](https://link.springer.com/article/10.1007/s42974-021-00042-x);
  - [*From zero to infinity: Minimum to maximum diversity of the planet
    by spatio-parametric Rao’s quadratic
    entropy*](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13270)
  - [*Helical graphs to visualize the NDVI temporal variation of forest
    vegetation in an open source
    space*](https://www.sciencedirect.com/science/article/abs/pii/S157495412200406X)
  - [*Phenological Diversity Trends with Remote Sensing
    Datacubes*](https://www.researchgate.net/profile/Matteo-Marcantonio/publication/381288639_Phenological_Diversity_Trends_with_Remote_Sensing_Datacubes/links/6666a0fbb769e76919269fa1/Phenological-Diversity-Trends-with-Remote-Sensing-Datacubes.pdf)
