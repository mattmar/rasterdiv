# Validate Input Parameters for Diversity Index Calculation

Validates the input parameters for diversity index calculation
functions. Checks for valid raster types, window sizes, alpha values,
and NA tolerance levels.

## Usage

``` r
validateInputs(x, window, alpha = 1, na.tolerance)
```

## Arguments

- x:

  Raster object to be validated.

- window:

  Size of the moving window for calculations.

- alpha:

  Diversity index parameter, default is 1.

- na.tolerance:

  Proportion of acceptable NA values within the window (range: 0 to 1).

## Value

None. Throws an error if any input is invalid.
