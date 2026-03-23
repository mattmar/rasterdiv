#' Parallelised Pielou's diversity index
#'
#' This function calculates Pielou's diversity index in a parallelized manner, 
#' allowing for improved performance on suitable hardware. The diversity index 
#' is computed using a moving window approach over the input data.
#'
#' @param x Input raster data, representing the environmental variable(s) 
#' over which the diversity index should be calculated.
#' @param window The size of the half-side of the square moving window used in the 
#' calculation. This determines the scale at which diversity is assessed.
#' @param na.tolerance A numeric value (between 0.0 and 1.0) indicating the 
#' proportion of NA values that are acceptable in each moving window over the 
#' raster data. If the proportion of NA values in a window exceeds this 
#' threshold, the resulting value for that window is set as NA. The default 
#' is 0.0, indicating no tolerance for NA values.
#' @param debugging Boolean flag indicating whether additional console
#' output should be generated for debugging purposes. Defaults to FALSE.
#' @param np The number of processes (cores) which will be spawned. Default value is 2.
#'
#' @return A matrix or list of matrices, depending on the input, containing 
#' the calculated Pielou diversity index values. Each cell in the output 
#' matrix represents the diversity index calculated from the corresponding 
#' moving window of the input data.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@gmail.com}, 
#' Martina Iannacito \email{martina.iannacito@inria.fr}, 
#' Duccio Rocchini \email{duccio.rocchini@unibo.it}
#'
#' @seealso \code{\link{Pielou}} for the non-parallelized version of the 
#' Pielou's diversity index calculation.
#'
#' @examples
#' \dontrun{
#' # Demonstration of function with hypothetical data
#' # Ensure you replace this with actual raster data
#' demo_raster <- #... (your raster data here)
#' result <- PielouP(x = demo_raster, win = 3, na.tolerance = 0.1, debugging = FALSE)
#' # proceed with analyzing 'result'
#' }

PielouP <- function(x, window = 1, na.tolerance = 1,
                    debugging = FALSE, np = 1, progBar = TRUE) {

  win <- window
  NAwin <- 2 * window + 1
  message("\n\nProcessing moving Window: ", NAwin)

  if (np > 1 && progBar) {
    message("Progress bar disabled for parallel execution.")
    progBar <- FALSE
  }

  values <- as.numeric(as.factor(x))
  x_1 <- matrix(data = values, nrow = nrow(x), ncol = ncol(x))

  hor <- matrix(NA, ncol = ncol(x), nrow = win)
  ver <- matrix(NA, ncol = win, nrow = nrow(x) + win * 2)
  tx <- cbind(ver, rbind(hor, x_1, hor), ver)

  rm(hor, ver, x_1, values)
  gc()

  PielouOP <- foreach::foreach(
    cl = (1 + win):(ncol(x) + win),
    .verbose = FALSE
  ) %dopar% {

    PielouOut <- sapply((1 + win):(nrow(x) + win), function(rw) {

      win_vals <- tx[(rw - win):(rw + win), (cl - win):(cl + win)]
      n_non_na <- sum(!is.na(win_vals))

      if (n_non_na < floor(NAwin^2 - ((NAwin^2) * na.tolerance))) {
        return(NA_real_)
      }

      tw <- summary(as.factor(win_vals), maxsum = 10000)
      if ("NA's" %in% names(tw)) {
        tw <- tw[-length(tw)]
      }

      if (debugging) {
        message(
          "\nPielou - parallelized\nWorking on coords ", rw, ",", cl,
          ". classes length: ", length(tw),
          ". window size = ", 2 * win + 1
        )
      }

      if (length(tw) <= 1) {
        return(0)
      }

      tw_values <- as.vector(tw)
      p <- tw_values / sum(tw_values)
      maxS <- log(length(tw))
      vv <- -sum(p * log(p)) / maxS

      if (debugging) {
        message(
          "\ncat: ", paste(names(tw), collapse = " "),
          " log S: ", maxS,
          " Pielou: ", vv
        )
      }

      return(vv)
    })

    PielouOut
  }

  message("\n\nParallel calculation of Pielou's index complete.\n")

  matrix(unlist(PielouOP), ncol = ncol(x), nrow = nrow(x), byrow = FALSE)
}