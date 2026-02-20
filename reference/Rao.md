# Rao's index

An alias for \`paRao\` with \`alpha\` fixed at 2.

## Usage

``` r
Rao(x, ...)
```

## Arguments

- x:

  Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster,
  or a list of these objects.

- ...:

  Other parameters passed to \`paRao\`.

## Value

A return value description.

## Examples

``` r
if (FALSE) { # \dontrun{
data(volcano)
r <- terra::rast(volcano)
res <- Rao(x = r, window = 3)
terra::plot(res[[1]][[1]])
} # }
```
