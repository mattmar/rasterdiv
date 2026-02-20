# Sequential Pielou's diversity index

Computes Pielou's diversity index using a sequential method,
particularly useful for handling large datasets that might not be
efficiently processed in a standard, non-sequential manner.

## Usage

``` r
PielouS(x, window = 1, na.tolerance = 1, debugging = FALSE)
```

## Arguments

- x:

  Input raster data, representing the environmental variable(s) over
  which the diversity index should be calculated.

- window:

  The size of the half-side of the square moving window used in the
  calculation. This determines the scale at which diversity is assessed.

- na.tolerance:

  A numeric value (between 0.0 and 1.0) indicating the proportion of NA
  values that are acceptable in each moving window over the raster data.
  If the proportion of NA values in a window exceeds this threshold, the
  resulting value for that window is set as NA. The default is 0.0,
  indicating no tolerance for NA values.

- debugging:

  Boolean flag indicating whether additional console output should be
  generated for debugging purposes. Defaults to FALSE.

## Value

A matrix or list of matrices, depending on the input, containing the
calculated Pielou diversity index values. Each cell in the output matrix
represents the diversity index calculated from the corresponding moving
window of the input data.

## See also

[`Pielou`](https://mattmar.github.io/rasterdiv/reference/Pielou.md) for
the standard computation of Pielou's diversity index.

## Author

Marcantonio Matteo <marcantoniomatteo@gmail.com>, Martina Iannacito
<martina.iannacito@inria.fr>, Duccio Rocchini <duccio.rocchini@unibo.it>

## Examples

``` r
if (FALSE) { # \dontrun{
# Demonstration of function with hypothetical data
# Ensure you replace this with actual raster data
demo_raster <- #... (your raster data here)
result <- PielouS(x = demo_raster, win = 3, na.tolerance = 0.1, debugging = FALSE)
# proceed with analyzing 'result'
} # }
```
