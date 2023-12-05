#' Hill's index of diversity - Hill numbers (D)
#'
#' Computes Hill's index of diversity (Hill numbers) on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects. In the latter case, only the first element of the list will be considered.
#' @param window The side of the square moving window. It must be an odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3.
#' @param alpha Order of the Hill number to compute the index. If \code{alpha} is a vector with length greater than 1, then the index will be calculated over \code{x} for each value in the sequence.
#' @param base The logarithm base for the calculation, default is natural logarithm.
#' @param rasterOut Boolean; if TRUE, the output will be in SpatRaster format with \code{x} as the template.
#' @param np The number of processes (cores) which will be spawned. Default value is 1.
#' @param na.tolerance A numeric value between 0.0 and 1.0, which indicates the proportion of NA values that will be tolerated to calculate Hill's index in each moving window over \code{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA; otherwise, Hill's index will be calculated considering the non-NA values. Default value is 1.0 (i.e., full tolerance for NA's).
#' @param cluster.type The type of cluster which will be created. Options are "MPI" (calls "makeMPIcluster"), "FORK," and "SOCK" (call "makeCluster"). Default type is "SOCK".
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed for debugging purposes.
#'
#' @details
#' Hill numbers (\eqn{{}^qD}) are calculated on numerical matrices as \eqn{{}^qD = (\sum_{i=1}^{R} {p^q}_i)^{1/(1-q)}}, where \emph{q} is the order of the Hill number, \emph{R} is the total number of categories (i.e., unique numerical values in a numerical matrix), and \emph{p} is the relative abundance of each category. When q=1, Shannon.R is called to calculate \eqn{exp(H^1)} instead of the indefinite \eqn{{}^1D}. If \eqn{q > 2*10^9}, BergerParker.R is called to calculate \eqn{1/{{}^\infty D}}. Hill numbers of low order weight more rare categories, whereas Hill numbers of higher order weight more dominant categories.
#'
#' @return A list of matrices of dimension \code{dim(x)} with length equal to the length of \code{alpha}.
#'
#' @references
#' Hill, M.O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology 54, 427-432.
#'
#' @note
#' Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
#' \code{apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi}
#' 
#' Microsoft Windows users may need some additional work to use "MPI". For more details, see:
#' \url{https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/}
#'
#' @seealso
#' \code{\link{BergerParker}}, \code{\link{Shannon}}
#'
#' @examples
#' # Minimal example; compute Hill's index with alpha 1:5 
#' a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
#' hill <- Hill(x=a,window=3,alpha=1:5)
#'
#' @export
Hill <- function(x, window = 3, alpha = 1, base = exp(1), rasterOut = TRUE, np = 1, na.tolerance=1, cluster.type = "SOCK", debugging = FALSE) {

  validateInputs(x, window, alpha, na.tolerance)
  rasterm <- prepareRaster(x)
  w <- calculateWindow(window)
  out <- if (np == 1) calculateHillSequential(rasterm[[1]], w, alpha, na.tolerance, debugging)
  else calculateHillParallel(rasterm[[1]], w, alpha, na.tolerance, debugging, base, cluster.type, np)
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
calculateHillSequential <- function(rasterm, w, alpha, na.tolerance, debugging) {
  if(debugging) {cat("#check: Before sequential function.")}
  lapply(w, function(win) {
    lapply(alpha, function(a) {
      outI <- if (abs(a - 1) < .Machine$double.eps) 
      exp(ShannonS(rasterm, win, na.tolerance, debugging))
      else if (a >= .Machine$integer.max) 
      1/BergerParkerS(rasterm, win, na.tolerance, debugging)
      else 
      HillS(rasterm, win, a, na.tolerance, debugging)
      })
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
calculateHillParallel <- function(rasterm, w, alpha, na.tolerance, debugging, base, cluster.type, np) {
  if(debugging) {cat("#check: Before parallel function.")}
  cls <- openCluster(cluster.type, np, debugging); on.exit(stopCluster(cls)); gc()
  lapply(w, function(win) {
    lapply(alpha, function(a) {
      outI <- if (abs(a - 1) < .Machine$double.eps) 
      exp(ShannonP(rasterm, window=win, na.tolerance=na.tolerance, debugging, np))
      else if (a >= .Machine$integer.max) 
      1/BergerParkerP(rasterm, window=win, na.tolerance=na.tolerance, debugging, np)
      else 
      HillP(rasterm, window=win, a, na.tolerance=na.tolerance, debugging, np)
      })
    })
}