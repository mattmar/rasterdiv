# Area-Based Sequential Parametric Rao's index of quadratic entropy (Q)

Calculates an area-based sequential version of the parametric Rao's
index of quadratic entropy (Q). This function is designed for situations
where the diversity index needs to consider geographical areas and works
with raster data representing the distribution of species or other
measures.

## Usage

``` r
mpaRaoAreaS(rasterm, area, alpha, simplify, dist_m, rescale, lambda, window)
```

## Arguments

- rasterm:

  Raster; the input raster data representing variables across a
  geographic space.

- area:

  Numeric; the input vector data representing the areas of interest.

- alpha:

  Numeric; alpha value for order of diversity in Hill's Index.

- simplify:

  Numeric; the parameter that determines the rounding off of the
  calculations.

- dist_m:

  Character; type of distance metric used (e.g., "euclidean",
  "manhattan", etc.).

- rescale:

  Logical; whether to scale and centre the values in each element of the
  raster data.

- lambda:

  Numeric; lambda parameter for Minkowski distance calculation.

- window:

  Numeric; defines the size of the moving window for the analysis.

## Value

A vector similar to the input, with additional columns representing
Rao's index values for each area.

## See also

[`paRao`](https://mattmar.github.io/rasterdiv/reference/paRao.md) for a
related function dealing with the parallel computation of Rao's index.

## Author

Matteo Marcantonio <marcantoniomatteo@gmail.com>, Duccio Rocchini
<duccio.rocchini@unibo.it>, Michele Torresani
<michele.torresani@unibo.it>
