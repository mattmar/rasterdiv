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

Pielou <- function(x, window = 3, rasterOut = TRUE, np = 1,
                   na.tolerance = 1, cluster.type = "SOCK",
                   debugging = FALSE, progBar = TRUE) {

  alpha <- 1
  validateInputs(x, window, alpha, na.tolerance)
  rasterm <- prepareRaster(x)
  w <- calculateWindow(window)

  out <- if (np == 1) {
    calculatePielouSequential(
      rasterm = rasterm[[1]],
      w = w,
      na.tolerance = na.tolerance,
      debugging = debugging,
      progBar = progBar
    )
  } else {
    calculatePielouParallel(
      rasterm = rasterm[[1]],
      w = w,
      na.tolerance = na.tolerance,
      debugging = debugging,
      cluster.type = cluster.type,
      np = np,
      progBar = progBar
    )
  }

  formatOutput(out, rasterOut, x, alpha, window)
}

calculatePielouSequential <- function(rasterm, w, na.tolerance, debugging, progBar = TRUE) {
  if (debugging) cat("#check: Before sequential function.")
  lapply(w, function(win) {
    PielouS(
      x = rasterm,
      window = win,
      na.tolerance = na.tolerance,
      debugging = debugging,
      progBar = progBar
    )
  })
}

calculatePielouParallel <- function(rasterm, w, na.tolerance, debugging,
                                    cluster.type, np, progBar = TRUE) {
  if (debugging) cat("#check: Before parallel function.")
  cls <- openCluster(cluster.type, np, debugging)
  on.exit(stopCluster(cls), add = TRUE)
  gc()

  lapply(w, function(win) {
    PielouP(
      x = rasterm,
      window = win,
      na.tolerance = na.tolerance,
      debugging = debugging,
      np = np,
      progBar = progBar
    )
  })
}
