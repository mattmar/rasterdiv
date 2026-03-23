#' Parallel Computation of Renyi's Diversity Index
#'
#' @description
#' This function computes Renyi's diversity index for each cell of a matrix, 
#' using a parallelized approach and considering a specified moving window.
#'
#' @param x A numeric matrix representing the data on which the index is to be calculated.
#' @param window The width of the moving window to consider for each cell. 
#'        The actual window size will be `(2 * window + 1) x (2 * window + 1)`. Default is 1.
#' @param alpha The alpha parameter for Renyi's index, influencing sensitivity 
#'        to species abundance. Default is 1.
#' @param base The base of the logarithm used in Renyi's formula. Default is `exp(1)` 
#'        (natural logarithm).
#' @param na.tolerance The tolerance level for missing data within the moving window. 
#'        A window will be processed only if the proportion of non-missing data is above this threshold. 
#'        Value should be between 0 and 1. Default is 1.
#' @param debugging Boolean flag to enable or disable debugging messages. Default is FALSE.
#' @param np Number of processes for parallel computation.#'
#' @return A matrix of the same dimensions as `x`, where each cell contains the 
#'         Renyi's diversity index calculated for the window around the cell.
#'
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' renyi_index <- RenyiP(data, window = 1, np = 1)
#'
#' @export

RenyiP <- function(x, window = 1, alpha = 1, base = exp(1),
                   na.tolerance = 1, debugging = FALSE, np = 1,
                   progBar = TRUE) {

  win <- window
  NAwin <- 2 * window + 1
  message("\n\nProcessing alpha: ", alpha, " Moving Window: ", NAwin)

  values <- as.numeric(as.factor(x))
  x_1 <- matrix(data = values, nrow = dim(x)[1], ncol = dim(x)[2])

  hor <- matrix(NA, ncol = dim(x)[2], nrow = win)
  ver <- matrix(NA, ncol = win, nrow = dim(x)[1] + win * 2)
  tx <- cbind(ver, rbind(hor, x_1, hor), ver)
  rm(hor, ver, x_1, values)
  gc()

  if (np > 1 && progBar) {
    message("Progress bar disabled for parallel execution.")
    progBar <- FALSE
  }

  RenyiOP <- foreach::foreach(
    cl = (1 + win):(dim(x)[2] + win),
    .verbose = FALSE
  ) %dopar% {

    if (debugging) {
      cat(cl)
    }

    RenyiOut <- sapply((1 + win):(dim(x)[1] + win), function(rw) {

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
          "Renyi - parallelized\nWorking on coords ", rw, ",", cl,
          ". classes length: ", length(tw),
          ". window size=", NAwin
        )
      }

      tw_values <- as.vector(tw)
      p <- tw_values / sum(tw_values)

      1 / (1 - alpha) * drop(log(sum(p^alpha), base))
    })

    RenyiOut
  }

  message("\n\nParallel calculation of Renyi's index complete.\n")

  matrix(unlist(RenyiOP), ncol = ncol(x), nrow = nrow(x), byrow = FALSE)
}