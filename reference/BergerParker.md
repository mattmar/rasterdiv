# Berger-Parker's diversity index

Computes Berger-Parker's diversity index on different classes of numeric
matrices using a moving window algorithm.

## Usage

``` r
BergerParker(
  x,
  window = 3,
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

  The side of the square moving window, it must be an odd numeric value
  greater than 1 to ensure that the target pixel is in the centre of the
  moving window. Default value is 3.

- rasterOut:

  Boolean, if TRUE, output will be in SpatRaster format with *x* as a
  template.

- np:

  The number of processes (cores) which will be spawned. Default value
  is 1.

- na.tolerance:

  A numeric value (0.0-1.0) which indicates the proportion of NA values
  that will be tolerated to calculate Berger-Parker's index in each
  moving window over *x*. If the relative proportion of NA's in a moving
  window is bigger than na.tolerance, then the value of the window will
  be set as NA, otherwise, Rao's index will be calculated considering
  the non-NA values. Default values are 1.0 (i.e., no tolerance for
  NA's).

- cluster.type:

  The type of cluster which will be created. The options are `"MPI"`
  (calls "makeMPIcluster"), `"FORK"`, and `"SOCK"` (call "makeCluster").
  Default type is `"SOCK"`.

- debugging:

  A boolean variable set to FALSE by default. If TRUE, additional
  messages will be printed. For de-bugging only.

## Value

A numerical matrix with dimensions as `dim(x)`.

## Details

Berger-Parker's index is the relative abundance of the most abundant
category (i.e., unique numerical values in the considered numerical
matrix). Berger-Parker's index equals the logarithm of the inverse
Renyi's index of order infinity, \\log(1/{}^\infty H)\\ or the inverse
of Hill's index of order infinity, \\1/{}^\infty D\\.

## Note

Linux users need to install libopenmpi for MPI parallel computing. Linux
Ubuntu users may try: apt-get update; apt-get upgrade; apt-get install
mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi

Microsoft Windows users may need some additional work to use "MPI",
see:  
<https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/>

## References

Berger, W.H., Parker, F.L. (1970). Diversity of planktonic foraminifera
in deep-sea sediments". Science, 168: 1345-1347.

## Author

Marcantonio Matteo <marcantoniomatteo@gmail.com>, Martina Iannacito
<martina.iannacito@inria.fr>, Duccio Rocchini <duccio.rocchini@unibo.it>

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal example; compute Renyi's index with alpha 1:5 
a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
berpar <- BergerParker(x=a, window=3)
} # }
```
