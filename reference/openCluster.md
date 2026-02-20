# Open a Parallel Cluster

Opens a parallel cluster for computation, registers it for parallel
operations, and ensures its closure on script exit.

## Usage

``` r
openCluster(cluster.type = "SOCK", np = 2, progBar = TRUE, debugging = FALSE)
```

## Arguments

- cluster.type:

  A character string specifying the type of cluster. Accepted values are
  "SOCK", "FORK", or "MPI".

- np:

  An integer specifying the number of processes to be used in the
  parallel cluster.

- progBar:

  logical. If TRUE a progress bar is shown.

- debugging:

  logical. For developer use.

## Value

An object representing the parallel cluster.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Open a SOCK cluster with 4 cores
  cls <- openCluster("SOCK", 4)
  # Your parallel computation code here
  # The cluster will automatically close when the script exits
} # }
```
