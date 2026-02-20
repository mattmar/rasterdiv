# Renyi Diversity Index Calculation

Computes Renyi diversity index for a given raster object. This function
allows specifying window size, alpha values, and various other
parameters for the calculation of the Renyi index.

## Usage

``` r
Renyi(
  x,
  window = 3,
  alpha = 1,
  base = exp(1),
  rasterOut = TRUE,
  np = 1,
  na.tolerance = 1,
  cluster.type = "SOCK",
  debugging = FALSE
)
```

## Arguments

- x:

  A raster object which can be a matrix, SpatialGridDataFrame,
  SpatRaster, list, or RasterStack.

- window:

  The size of the moving window; must be an odd integer.

- alpha:

  A numeric vector of alpha values for the Renyi index.

- base:

  The logarithm base for the calculation, default is natural logarithm.

- rasterOut:

  Logical; if TRUE, returns a SpatRaster object, otherwise returns a
  list.

- np:

  Number of processes for parallel computation.

- na.tolerance:

  Tolerance level for NA values, must be within \[0-1\].

- cluster.type:

  Type of cluster for parallel computation, either "SOCK" or "MPI".

- debugging:

  Logical; if TRUE, provides additional console output for debugging.

## Value

A SpatRaster object or a list of calculated Renyi indices.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- Renyi(ndvi.8bit, window = 3, alpha = c(0, 1, 2))
} # }
```
