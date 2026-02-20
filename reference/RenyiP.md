# Parallel Computation of Renyi's Diversity Index

This function computes Renyi's diversity index for each cell of a
matrix, using a parallelized approach and considering a specified moving
window.

## Usage

``` r
RenyiP(
  x,
  window = 1,
  alpha = 1,
  base = exp(1),
  na.tolerance = 1,
  debugging = FALSE,
  np = 1
)
```

## Arguments

- x:

  A numeric matrix representing the data on which the index is to be
  calculated.

- window:

  The width of the moving window to consider for each cell. The actual
  window size will be \`(2 \* window + 1) x (2 \* window + 1)\`. Default
  is 1.

- alpha:

  The alpha parameter for Renyi's index, influencing sensitivity to
  species abundance. Default is 1.

- base:

  The base of the logarithm used in Renyi's formula. Default is
  \`exp(1)\` (natural logarithm).

- na.tolerance:

  The tolerance level for missing data within the moving window. A
  window will be processed only if the proportion of non-missing data is
  above this threshold. Value should be between 0 and 1. Default is 1.

- debugging:

  Boolean flag to enable or disable debugging messages. Default is
  FALSE.

- np:

  Number of processes for parallel computation.#'

## Value

A matrix of the same dimensions as \`x\`, where each cell contains the
Renyi's diversity index calculated for the window around the cell.

## Examples

``` r
data <- matrix(runif(100), nrow = 10)
renyi_index <- RenyiP(data, window = 1, np = 1)
#> 
#> 
#> Processing alpha: 1 Moving Window: 3
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> 
#> 
#>  Parallel calculation of Renyi's index complete.
```
