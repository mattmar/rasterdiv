# Hill's index of diversity - Hill numbers (D)

Computes Hill's index of diversity (Hill numbers) on different classes
of numeric matrices using a moving window algorithm.

## Usage

``` r
Hill(
  x,
  window = 3,
  alpha = 1,
  base = exp(1),
  rasterOut = TRUE,
  np = 1,
  na.tolerance = 1,
  cluster.type = "SOCK",
  debugging = FALSE
)
```

## Arguments

- x:

  Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster,
  or a list of these objects. In the latter case, only the first element
  of the list will be considered.

- window:

  The side of the square moving window. It must be an odd numeric value
  greater than 1 to ensure that the target pixel is in the centre of the
  moving window. Default value is 3.

- alpha:

  Order of the Hill number to compute the index. If `alpha` is a vector
  with length greater than 1, then the index will be calculated over `x`
  for each value in the sequence.

- base:

  The logarithm base for the calculation, default is natural logarithm.

- rasterOut:

  Boolean; if TRUE, the output will be in SpatRaster format with `x` as
  the template.

- np:

  The number of processes (cores) which will be spawned. Default value
  is 1.

- na.tolerance:

  A numeric value between 0.0 and 1.0, which indicates the proportion of
  NA values that will be tolerated to calculate Hill's index in each
  moving window over `x`. If the relative proportion of NA's in a moving
  window is bigger than na.tolerance, then the value of the window will
  be set as NA; otherwise, Hill's index will be calculated considering
  the non-NA values. Default value is 1.0 (i.e., full tolerance for
  NA's).

- cluster.type:

  The type of cluster which will be created. Options are "MPI" (calls
  "makeMPIcluster"), "FORK," and "SOCK" (call "makeCluster"). Default
  type is "SOCK".

- debugging:

  A boolean variable set to FALSE by default. If TRUE, additional
  messages will be printed for debugging purposes.

## Value

A list of matrices of dimension `dim(x)` with length equal to the length
of `alpha`.

## Details

Hill numbers (\\{}^qD\\) are calculated on numerical matrices as \\{}^qD
= (\sum\_{i=1}^{R} {p^q}\_i)^{1/(1-q)}\\, where *q* is the order of the
Hill number, *R* is the total number of categories (i.e., unique
numerical values in a numerical matrix), and *p* is the relative
abundance of each category. When q=1, Shannon.R is called to calculate
\\exp(H^1)\\ instead of the indefinite \\{}^1D\\. If \\q \> 2\*10^9\\,
BergerParker.R is called to calculate \\1/{{}^\infty D}\\. Hill numbers
of low order weight more rare categories, whereas Hill numbers of higher
order weight more dominant categories.

## Note

Linux users need to install libopenmpi for MPI parallel computing. Linux
Ubuntu users may try:
`apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi`

Microsoft Windows users may need some additional work to use "MPI". For
more details, see:
<https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/>

## References

Hill, M.O. (1973). Diversity and evenness: a unifying notation and its
consequences. Ecology 54, 427-432.

## See also

[`BergerParker`](https://mattmar.github.io/rasterdiv/reference/BergerParker.md),
[`Shannon`](https://mattmar.github.io/rasterdiv/reference/Shannon.md)

## Examples

``` r
# Minimal example; compute Hill's index with alpha 1:5 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
hill <- Hill(x=a,window=3,alpha=1:5)
#> 
#> 
#> Processing moving Window: 3
#> 
#> 
#> Processing alpha: 2 Moving Window: 3
#> 
#> 
#> Processing alpha: 3 Moving Window: 3
#> 
#> 
#> Processing alpha: 4 Moving Window: 3
#> 
#> 
#> Processing alpha: 5 Moving Window: 3
```
