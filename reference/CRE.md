# Cumulative Residual Entropy (CRE) Function

Computes the Cumulative Residual Entropy (CRE) for spatial raster data.
This function can be used with either a single raster layer or a list of
raster layers. It supports both classic and multidimensional methods for
CRE computation.

## Usage

``` r
CRE(
  x,
  window = 3,
  method = "classic",
  rasterOut = TRUE,
  rescale = FALSE,
  na.tolerance = 1,
  simplify = 2,
  np = 1,
  cluster.type = "SOCK",
  progBar = TRUE,
  debugging = FALSE
)
```

## Arguments

- x:

  A matrix, SpatRaster, or a list of SpatRaster objects.

- window:

  The size of the moving window, must be an odd integer.

- method:

  The method for CRE computation, either "classic" or
  "multidimensional".

- rasterOut:

  Logical, if TRUE, returns a SpatRaster, else returns a matrix.

- rescale:

  Logical, if TRUE, rescales the data before processing.

- na.tolerance:

  A numeric value between 0 and 1, indicating the tolerance level for NA
  values.

- simplify:

  Integer, the number of decimal places for data rounding in case of
  float numbers.

- np:

  The number of parallel processes to use.

- cluster.type:

  The type of parallel cluster to use, options are "SOCK", "FORK", or
  "MPI".

- progBar:

  logical. If TRUE a progress bar is shown.

- debugging:

  Logical, if TRUE, provides additional debugging information during
  execution.

## Value

Depending on the 'rasterOut' parameter, this function returns either a
SpatRaster or a matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
# For a matrix input:
result <- CRE(matrix_data, window=3, method="classic")

# For a SpatRaster input:
result <- CRE(raster_data, window=3, method="classic", rasterOut=TRUE)
} # }
```
