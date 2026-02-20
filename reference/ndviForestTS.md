# Simulated NDVI dataset

A `list` of 8-bit matrices.

## Format

A `list` containing matrices:

- ndviForestTS:

  List of matrixes of 9 cells simulating NDVI of a patch of forests over
  3 years. Each matrix represents a day in the time series.

## Details

This list represents a time series of NDVI values of a patch of forest
over 3 years. It is stored as a `list`, suitable for explaining how to
make helical plots.

## Examples

``` r
ndviForestTS <- readRDS(system.file("extdata", "ndviForestTS.rds", package = "rasterdiv"))
```
