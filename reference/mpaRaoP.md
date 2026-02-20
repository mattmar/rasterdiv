# Multidimensional parallel Parametric Rao's index of quadratic entropy (Q)

Multidimensional parametric Rao's index of quadratic entropy (Q).

## Usage

``` r
mpaRaoP(
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

  alpha value for order of diversity in Hill's Index.

- window:

  half of the side of the square moving window.

- dist_m:

  Type of distance used.

- na.tolerance:

  a numeric value \\(0.0-1.0)\\ which indicates the proportion of NA
  values that will be tolerated to calculate Rao's index in each moving
  window over *x*. If the relative proportion of NA's in a moving window
  is bigger than na.tolerance, then the value of the window will be set
  as NA, otherwise Rao's index will be calculated considering the non-NA
  values. Default values is 0.0 (i.e., no tolerance for NA's).

- rescale:

  Scale and centre values in each of the element of x.

- lambda:

  Lambda value for Minkowski distance.

- diag:

  Boolean. Diagonal of the distance matrix.

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

  a boolean variable set to FALSE by default. If TRUE, additional
  messages will be printed. For de-bugging only.

- isfloat:

  Are the input data floats?

- mfactor:

  Multiplication factor in case of input data as float numbers.

- np:

  the number of processes (cores) which will be spawned.

- progBar:

  logical. If TRUE a progress bar is shown.

## Value

A list of matrices of dimension `dim(x)` with length equal to the length
of `alpha`.

## See also

[`paRao`](https://mattmar.github.io/rasterdiv/reference/paRao.md)

## Author

Duccio Rocchini <duccio.rocchini@unibo.it>, Marcantonio Matteo
<marcantoniomatteo@gmail.com>
