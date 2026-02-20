# Parallelised Hill's diversity index

Parallelised computation of Hill's diversity index.

## Usage

``` r
HillP(rasterm, window, alpha, na.tolerance, debugging, np)
```

## Arguments

- rasterm:

  Input data.

- window:

  Half of the side of the square moving window.

- alpha:

  Alpha value for the order of diversity in Hill's Index.

- na.tolerance:

  A numeric value between 0.0 and 1.0, which indicates the proportion of
  NA values that will be tolerated to calculate the index in each moving
  window over `rasterm`. If the relative proportion of NA's in a moving
  window is bigger than na.tolerance, then the value of the window will
  be set as NA; otherwise, the index will be calculated considering the
  non-NA values. Default value is 0.0 (i.e., no tolerance for NA's).

- debugging:

  A boolean variable set to FALSE by default. If TRUE, additional
  messages will be printed for debugging purposes.

- np:

  Number of processes for parallel computation.

## Value

Matrix or a list of matrices with the Hill index computed through a
moving window of the given size.

## See also

[`Hill`](https://mattmar.github.io/rasterdiv/reference/Hill.md)

## Author

Marcantonio Matteo <marcantoniomatteo@gmail.com>, Martina Iannacito
<martina.iannacito@inria.fr>, Duccio Rocchini <duccio.rocchini@unibo.it>
