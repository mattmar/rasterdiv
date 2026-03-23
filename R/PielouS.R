#' Sequential Pielou's diversity index
#'
#' Computes Pielou's diversity index using a sequential method, particularly useful 
#' for handling large datasets that might not be efficiently processed in a 
#' standard, non-sequential manner.
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
#' @seealso \code{\link{Pielou}} for the standard computation of Pielou's 
#' diversity index.
#'
#' @examples
#' \dontrun{
#' # Demonstration of function with hypothetical data
#' # Ensure you replace this with actual raster data
#' demo_raster <- #... (your raster data here)
#' result <- PielouS(x = demo_raster, win = 3, na.tolerance = 0.1, debugging = FALSE)
#' # proceed with analyzing 'result'
#' }

PielouS <- function(x, window = 1, na.tolerance = 1, debugging = FALSE, progBar = TRUE) {

  win <- window
  NAwin <- 2 * window + 1
  message("\n\nProcessing moving Window: ", NAwin)

  if (progBar) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent in :elapsed",
      total = ncol(x),
      clear = FALSE,
      width = 60,
      force = FALSE
    )
  }

  out <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x))

  values <- as.numeric(as.factor(x))
  x_1 <- matrix(data = values, nrow = nrow(x), ncol = ncol(x))

  # Add additional columns and rows for moving window
  hor <- matrix(NA, ncol = ncol(x), nrow = win)
  ver <- matrix(NA, ncol = win, nrow = nrow(x) + win * 2)
  tx <- cbind(ver, rbind(hor, x_1, hor), ver)

  for (cl in (1 + win):(ncol(x) + win)) {
    if (progBar) pb$tick()

    for (rw in (1 + win):(nrow(x) + win)) {

      win_vals <- tx[(rw - win):(rw + win), (cl - win):(cl + win)]
      n_non_na <- sum(!is.na(win_vals))

      if (n_non_na < floor(NAwin^2 - ((NAwin^2) * na.tolerance))) {
        out[rw - win, cl - win] <- NA_real_
      } else {
        tw <- summary(as.factor(win_vals), maxsum = 10000)

        if ("NA's" %in% names(tw)) {
          tw <- tw[-length(tw)]
        }

        if (debugging) {
          message(
            "\nPielou\nWorking on coords ", rw, ",", cl,
            ". classes length: ", length(tw),
            ". window size = ", NAwin
          )
        }

        if (length(tw) <= 1) {
          out[rw - win, cl - win] <- 0
        } else {
          tw_values <- as.vector(tw)
          p <- tw_values / sum(tw_values)
          maxS <- log(length(tw))
          out[rw - win, cl - win] <- -sum(p * log(p)) / maxS
        }

        if (debugging) {
          message(
            "\ncat: ", paste(names(tw), collapse = " "),
            " log S: ", if (length(tw) > 1) log(length(tw)) else 0,
            " Pielou: ", out[rw - win, cl - win]
          )
        }
      }
    }
  }

  return(out)
}