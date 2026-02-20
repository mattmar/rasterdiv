# Sequential Berger-Parker's diversity index

This function calculates the Berger-Parker's diversity index for each
cell in a matrix, considering a specified moving window around each
cell.

## Usage

``` r
BergerParkerS(x, window = 1, na.tolerance = 1, debugging = FALSE)
```

## Arguments

- x:

  A numeric matrix representing the data on which the index is to be
  calculated.

- window:

  The width of the moving window to consider for each cell. The actual
  window size will be \`(2 \* window + 1) x (2 \* window + 1)\`. Default
  is 1.

- na.tolerance:

  The tolerance level for missing data within the moving window. A
  window will be processed only if the proportion of non-missing data is
  above this threshold. Value should be between 0 and 1. Default is 1.

- debugging:

  Boolean flag to enable or disable debugging messages. Default is
  FALSE.

## Value

A matrix of the same dimensions as \`x\`, where each cell contains the
Berger-Parker's diversity index calculated for the window around the
cell.

## Details

Berger-Parker's diversity index calculated sequentially over a raster
matrix.

## Examples

``` r
data <- matrix(runif(100), nrow = 10)
bp_index <- BergerParkerS(data, window = 1)
#> 
#> 
#> Processing moving Window: 3
```
