#' Calculate Berger-Parker Index on a Matrix
#'
#' @description
#' This function computes Berger-Parker Index for each cell of a matrix, 
#' using a parallelized approach and considering a specified moving window.
#'
#' @param x A numeric matrix representing the data on which the index is to be calculated.
#' @param window The width of the moving window to consider for each cell. 
#'        The actual window size will be `(2 * window + 1) x (2 * window + 1)`. Default is 1.
#' @param na.tolerance The tolerance level for missing data within the moving window. 
#'        A window will be processed only if the proportion of non-missing data is above this threshold. 
#'        Value should be between 0 and 1. Default is 1.
#' @param debugging Boolean flag to enable or disable debugging messages. Default is FALSE.
#' @param np The number of processes (cores) which will be spawned. Default value is 2.
#'
#' @return A matrix of the same dimensions as `x`, where each cell contains the 
#'         Berger-Parker Index calculated for the window around the cell.
#'
#' @examples
#' \dontrun{
#' data <- matrix(runif(100), nrow = 10)
#' bp_index <- BergerParkerP(data, window = 1, np=2)
#' }
#'
#' @export
BergerParkerP <- function(x, window = 1, na.tolerance = 1,
                          debugging = FALSE, np = 1, progBar = TRUE) {

  win <- window
  NAwin <- 2 * window + 1
  message("\n\nProcessing moving Window: ", NAwin)

  if (np > 1 && progBar) {
    message("Progress bar disabled for parallel execution.")
    progBar <- FALSE
  }

  # Reshape values
  values <- as.numeric(as.factor(x))
  x_1 <- matrix(values, nrow = nrow(x), ncol = ncol(x))

  # Add padding to match moving window
  hor <- matrix(NA, ncol = ncol(x), nrow = win)
  ver <- matrix(NA, ncol = win, nrow = nrow(x) + win * 2)
  tx <- cbind(ver, rbind(hor, x_1, hor), ver)

  rm(hor, ver, x_1, values)
  gc()

  BergerParkerOP <- foreach::foreach(
    cl = (1 + win):(ncol(x) + win),
    .verbose = FALSE
  ) %dopar% {

    if (debugging) {
      cat(cl)
    }

    BergerParkerOut <- sapply((1 + win):(nrow(x) + win), function(rw) {

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
          "Berger-Parker - parallelised\nWorking on coords ", rw, ",", cl,
          ". classes length: ", length(tw),
          ". window size=", NAwin
        )
      }

      tw_values <- as.vector(tw)
      max(tw_values / sum(tw_values))
    })

    BergerParkerOut
  }

  message("\n\nParallel calculation of Berger Parker's index complete.\n")

  matrix(unlist(BergerParkerOP), ncol = ncol(x), nrow = nrow(x), byrow = FALSE)
}