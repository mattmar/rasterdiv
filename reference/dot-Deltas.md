# Calculate Differences Among Values

This function computes the differences among values of a table, used in
probability calculations.

## Usage

``` r
.Deltas(P, first = 0)
```

## Arguments

- P:

  A numeric vector or matrix representing probabilities.

- first:

  The starting value for difference calculation.

## Value

A vector or matrix of differences.

## Examples

``` r
P <- c(0.2, 0.3, 0.5)
.Deltas(P)
#> [1] 1
```
