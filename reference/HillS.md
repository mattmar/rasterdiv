# Sequential Hill's diversity index

Computes Hill's diversity index in a non-parallelized (sequential)
manner.

## Usage

``` r
HillS(x, window = 1, alpha = 1, na.tolerance = 1, debugging = FALSE)
```

## Arguments

- x:

  Input data.

- window:

  Half of the side of the square moving window.

- alpha:

  Alpha value for the order of diversity in Hill's Index.

- na.tolerance:

  A numeric value between 0.0 and 1.0, which indicates the proportion of
  NA values that will be tolerated to calculate Hill's index in each
  moving window over `x`. If the relative proportion of NA's in a moving
  window is greater than na.tolerance, then the value of the window will
  be set as NA; otherwise, the index will be calculated considering the
  non-NA values. Default value is 1.0 (i.e., full tolerance for NA's).

- debugging:

  A boolean variable set to FALSE by default. If TRUE, additional
  messages will be printed for debugging purposes.

## Value

Matrix or a list of matrices with the Hill index computed through a
moving window of the given size.

## See also

[`Hill`](https://mattmar.github.io/rasterdiv/reference/Hill.md)

## Author

Marcantonio Matteo <marcantoniomatteo@gmail.com>, Martina Iannacito
<martina.iannacito@inria.fr>, Duccio Rocchini <duccio.rocchini@unibo.it>
