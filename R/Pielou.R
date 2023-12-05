#' Pielou's Evenness Index
#'
#' @description Calculates Pielou's Evenness Index for a given raster object over a specified window size. The function can operate in either sequential or parallel mode.
#'
#' @param x A raster object (matrix, SpatRaster, SpatialGridDataFrame, or a list containing one of these).
#' @param window The size of the moving window to be used for the calculation. Must be an odd integer.
#' @param rasterOut Logical, if TRUE the output will be a raster object; if FALSE a matrix.
#' @param np The number of processes to use in parallel mode. If np > 1, parallel computation is enabled.
#' @param na.tolerance The tolerance level for NA values within the moving window, expressed as a proportion (0 to 1).
#' @param cluster.type The type of cluster to use for parallel computation (e.g., "SOCK", "FORK").
#' @param debugging Logical, if TRUE debugging information will be printed.
#'
#' @return Returns a raster object or matrix containing the Pielou's Evenness Index values.
#' @export

Pielou <- function(x, window = 3, rasterOut = TRUE, np = 1, na.tolerance=1, cluster.type = "SOCK", debugging = FALSE) {

  alpha=1 
  validateInputs(x, window, alpha, na.tolerance)
  rasterm <- prepareRaster(x)
  w <- calculateWindow(window)
  out <- if (np == 1) calculatePielouSequential(rasterm[[1]], w, na.tolerance, debugging)
  else calculatePielouParallel(rasterm[[1]], w, na.tolerance, debugging, cluster.type, np)
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

calculatePielouSequential <- function(rasterm, w, na.tolerance, debugging) {
  if(debugging) {cat("#check: Before sequential function.")}
  lapply(w, function(win) {
    PielouS(rasterm, win, na.tolerance, debugging)
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
calculatePielouParallel <- function(rasterm, w, na.tolerance, debugging, cluster.type, np) {
  if(debugging) {cat("#check: Before parallel function.")}
  cls <- openCluster(cluster.type, np, debugging); on.exit(stopCluster(cls)); gc()
  lapply(w, function(win) {
    PielouP(rasterm, win, na.tolerance, debugging, np)
    })
}
