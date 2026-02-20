# Cumulative Residual Entropy

This function calculates the Cumulative Residual Entropy (CRE) for a
given set of values.

## Usage

``` r
.CRE_(B, base = exp(1))
```

## Arguments

- B:

  A numeric vector or matrix representing the values for which CRE is to
  be calculated.

- base:

  The base of the logarithm used in the calculation. The default is the
  natural logarithm (e).

## Value

A numeric value representing the CRE.

## Examples

``` r
B <- c(1, 2, 3, 4)
.CRE_(B)
#> [1] 0.9089087
```
