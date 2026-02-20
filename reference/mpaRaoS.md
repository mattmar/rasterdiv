# Multidimensional sequential Parametric Rao's index of quadratic entropy (Q)

This function calculates the multidimensional parametric Rao's index of
quadratic entropy (Q) using a sequential method. It is particularly
useful in contexts where parallel computation is not feasible or
desired. The function applies a moving window approach to the provided
raster data stack.

## Usage

``` r
mpaRaoS(
  x,
  alpha,
  window,
  dist_m,
  na.tolerance,
  rescale,
  lambda,
  diag,
  time_vector,
  stepness,
  midpoint,
  cycle_length,
  time_scale,
  debugging,
  isfloat,
  mfactor,
  np,
  progBar
)
```

## Arguments

- x:

  input list.

- alpha:

  Numeric; alpha value for order of diversity in Hill's Index.

- window:

  Numeric; half of the side of the square moving window used for
  calculation.

- dist_m:

  Character; type of distance used in the analysis.

- na.tolerance:

  Numeric; a threshold between 0.0 and 1.0 indicating the allowable
  proportion of NA values within each moving window. If the proportion
  of NA values exceeds this, the window's value is set as NA; otherwise,
  the computation uses the non-NA values.

- rescale:

  Logical; if TRUE, scales and centres the values in each element of
  'x'.

- lambda:

  Numeric; lambda value used for Minkowski distance calculation.

- diag:

  Logical; if TRUE, includes the diagonal of the distance matrix in
  computations.

- time_vector:

  time;

- stepness:

  numeric; steepness of the logistic function.

- midpoint:

  numeric; midpoint of the logistic function

- cycle_length:

  string; The length of the cycle. Can be a numeric value or a string
  specifying the units ('year', 'month', 'day', 'hour', 'minute',
  'second'). When numeric, the cycle length is in the same units as
  time_scale. When a string, it specifies the time unit of the cycle.

- time_scale:

  string; Specifies the time scale for the conversion. Must be one of
  'year', 'month', 'day', 'hour', 'minute', 'second'. When cycle_length
  is a string, time_scale changes the unit in which the result is
  expressed. When cycle_length is numeric, time_scale is used to compute
  the elapsed time in seconds.

- debugging:

  Logical; if TRUE, additional diagnostic messages are output, useful
  for debugging. Default is FALSE.

- isfloat:

  Logical; specifies if the input data are floats.

- mfactor:

  Numeric; multiplication factor applied if input data are float
  numbers.

- np:

  Number of processes for parallel computation.

- progBar:

  logical. If TRUE a progress bar is shown.

## Value

A list of matrices, each representing a layer of the input RasterStack,
containing calculated Rao's index values. The dimensions correspond to
those of the input, and the list length is equal to the length of
'alpha'.

## See also

[`paRao`](https://mattmar.github.io/rasterdiv/reference/paRao.md) for
the parallelized version of the Rao's index computation.

## Author

Duccio Rocchini <duccio.rocchini@unibo.it>, Matteo Marcantonio
<marcantoniomatteo@gmail.com>
