# Calculate Cumulative Residual Probability

This function computes the cumulative residual probability for a given
set of probabilities.

## Usage

``` r
.CumRes(a)
```

## Arguments

- a:

  A numeric vector or matrix representing probabilities.

## Value

A numeric vector or matrix of cumulative residual probabilities.

## Examples

``` r
a <- data.frame(V1= c(0.2, 0.3, 0.5), V2 =c(0.2, 0.3, 0.5))
.CumRes(a)
#>       V1  V2
#> [1,] 2.0 1.0
#> [2,] 1.6 0.8
#> [3,] 1.0 0.5
```
