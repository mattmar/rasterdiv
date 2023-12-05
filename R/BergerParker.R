#' Berger-Parker's diversity index
#'
#' Computes Berger-Parker's diversity index on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects. In the latter case, only the first element of the list will be considered.
#' @param window The side of the square moving window, it must be an odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3.
#' @param rasterOut Boolean, if TRUE, output will be in SpatRaster format with \emph{x} as a template.
#' @param np The number of processes (cores) which will be spawned. Default value is 1.
#' @param na.tolerance A numeric value (0.0-1.0) which indicates the proportion of NA values that will be tolerated to calculate Berger-Parker's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise, Rao's index will be calculated considering the non-NA values. Default values are 1.0 (i.e., no tolerance for NA's).
#' @param cluster.type The type of cluster which will be created. The options are \code{"MPI"} (calls "makeMPIcluster"), \code{"FORK"}, and \code{"SOCK"} (call "makeCluster"). Default type is \code{"SOCK"}.
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For de-bugging only.
#' @details Berger-Parker's index is the relative abundance of the most abundant category (i.e., unique numerical values in the considered numerical matrix). Berger-Parker's index equals the logarithm of the inverse Renyi's index of order infinity, \eqn{log(1/{}^\infty H)} or the inverse of Hill's index of order infinity, \eqn{1/{}^\infty D}.
#' @return A numerical matrix with dimensions as \code{dim(x)}.
#' @references Berger, W.H., Parker, F.L. (1970). Diversity of planktonic foraminifera in deep-sea sediments". Science, 168: 1345-1347.
#' @author Marcantonio Matteo \email{marcantoniomatteo@gmail.com}, Martina Iannacito \email{martina.iannacito@inria.fr}, Duccio Rocchini \email{duccio.rocchini@unibo.it}
#' @note Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
#' apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi
#'
#' Microsoft Windows users may need some additional work to use "MPI", see:\cr
#' \url{https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/}
#'
#' @examples
#' \dontrun{
#' # Minimal example; compute Renyi's index with alpha 1:5 
#' a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
#' berpar <- BergerParker(x=a, window=3)
#' }
#' @export
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach %dopar%

BergerParker <- function(x, window = 3, rasterOut = TRUE, np = 1, na.tolerance=1, cluster.type = "SOCK", debugging = FALSE) {

  alpha=1 
  validateInputs(x, window, alpha, na.tolerance)
  rasterm <- prepareRaster(x)
  w <- calculateWindow(window)
  out <- if (np == 1) calculateBergerParkerSequential(rasterm[[1]], w, na.tolerance, debugging)
  else calculateBergerParkerParallel(rasterm[[1]], w, na.tolerance, debugging, cluster.type, np)
  formatOutput(out, rasterOut, x, alpha, window)
}

#' Calculate Sequentially
#'
#' @description Internal function to calculate indices sequentially.
#'
#' @param rasterm Prepared raster object for computation.
#' @param w The operative moving window size.
#' @param alpha The alpha parameter (used in some indices).
#' @param na.tolerance Proportion of acceptable NA values in the window.
#' @param debugging Logical flag for debugging mode.
#'
#' @return Returns a list or matrix of calculated index values.
#' @noRd

calculateBergerParkerSequential <- function(rasterm, w, na.tolerance, debugging) {
  if(debugging) {cat("#check: Before sequential function.")}
  lapply(w, function(win) {
    BergerParkerS(rasterm, win, na.tolerance, debugging)
    })
}

#' Calculate in Parallel
#'
#' @description Internal function to calculate indices in parallel.
#'
#' @param rasterm Prepared raster object for computation.
#' @param w The operative moving window size.
#' @param alpha The alpha parameter (used in some indices).
#' @param na.tolerance Proportion of acceptable NA values in the window.
#' @param debugging Logical flag for debugging mode.
#' @param cluster.type Cluster type for parallel computation.
#' @param np Number of processes for parallel computation.
#'
#' @return Returns a list or matrix of calculated index values.
#' @noRd
calculateBergerParkerParallel <- function(rasterm, w, na.tolerance, debugging, cluster.type, np) {
  if(debugging) {cat("#check: Before parallel function.")}
    cls <- openCluster(cluster.type, np, debugging); on.exit(stopCluster(cls)); gc()
  lapply(w, function(win) {
    BergerParkerP(rasterm, win, na.tolerance, debugging, np)
    })
}
