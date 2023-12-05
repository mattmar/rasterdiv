#' Renyi Diversity Index Calculation
#'
#' @description
#' Computes Renyi diversity index for a given raster object. This function allows
#' specifying window size, alpha values, and various other parameters for 
#' the calculation of the Renyi index.
#'
#' @param x A raster object which can be a matrix, SpatialGridDataFrame, SpatRaster, list, or RasterStack.
#' @param window The size of the moving window; must be an odd integer.
#' @param alpha A numeric vector of alpha values for the Renyi index.
#' @param base The logarithm base for the calculation, default is natural logarithm.
#' @param rasterOut Logical; if TRUE, returns a SpatRaster object, otherwise returns a list.
#' @param np Number of processes for parallel computation.
#' @param na.tolerance Tolerance level for NA values, must be within [0-1].
#' @param cluster.type Type of cluster for parallel computation, either "SOCK" or "MPI".
#' @param debugging Logical; if TRUE, provides additional console output for debugging.
#'
#' @return A SpatRaster object or a list of calculated Renyi indices.
#'
#' @examples
#' \dontrun{
# Assuming ndvi.8bit is a list of SpatRaster objects
#' result <- Renyi(ndvi.8bit, window = 3, alpha = c(0, 1, 2))
#' }
#'
#' @export
Renyi <- function(x, window = 3, alpha = 1, base = exp(1), rasterOut = TRUE, np = 1, na.tolerance=1, cluster.type = "SOCK", debugging = FALSE) {

  validateInputs(x, window, alpha, na.tolerance)
  rasterm <- prepareRaster(x)
  w <- calculateWindow(window)
  out <- if (np == 1) calculateRenyiSequential(rasterm[[1]], w, alpha, na.tolerance, debugging, base)
  else calculateRenyiParallel(rasterm[[1]], w, alpha, na.tolerance, debugging, base, cluster.type, np)
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
calculateRenyiSequential <- function(rasterm, w, alpha, na.tolerance, debugging, base) {
  if(debugging) {cat("#check: Before sequential function.")}
  lapply(w, function(win) {
    lapply(alpha, function(a) {
      outI <- if (abs(a - 1) < .Machine$double.eps) 
      ShannonS(rasterm, win, na.tolerance, debugging)
      else if (a >= .Machine$integer.max) 
      BergerParkerS(rasterm, win, na.tolerance, debugging)
      else 
      RenyiS(rasterm, win, a, base, na.tolerance, debugging)
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
calculateRenyiParallel <- function(rasterm, w, alpha, na.tolerance, debugging, base, cluster.type, np) {
  if(debugging) {cat("#check: Before parallel function.")}
  cls <- openCluster(cluster.type, np, debugging); on.exit(stopCluster(cls)); gc()
  lapply(w, function(win) {
    lapply(alpha, function(a) {
      outI <- if (abs(a - 1) < .Machine$double.eps) 
      ShannonP(rasterm, win, na.tolerance, debugging, np)
      else if (a >= .Machine$integer.max) 
      BergerParkerP(rasterm, win, na.tolerance, debugging, np)
      else 
      RenyiP(rasterm, win, a, base, na.tolerance, debugging, np)
      })
    })
}