# Calculate Point Probability

This function computes the probability of each point in a given vector
or matrix.

## Usage

``` r
.Prob(C)
```

## Arguments

- C:

  A numeric vector or matrix.

## Value

A vector of probabilities corresponding to each point in \`C\`.

## Examples

``` r
C <- c(1, 1, 2, 2, 3)
.Prob(C)
#> C
#>   1   2   3 
#> 0.4 0.4 0.2 
```
