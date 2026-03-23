#' Multidimensional sequential Parametric Rao's index of quadratic entropy (Q)
#'
#' This function calculates the multidimensional parametric Rao's index of quadratic entropy (Q) using a 
#' sequential method. It is particularly useful in contexts where parallel computation is not feasible or desired.
#' The function applies a moving window approach to the provided raster data stack.
#'
#' @param x input list.
#' @param alpha Numeric; alpha value for order of diversity in Hill's Index.
#' @param window Numeric; half of the side of the square moving window used for calculation.
#' @param dist_m Character; type of distance used in the analysis.
#' @param na.tolerance Numeric; a threshold between 0.0 and 1.0 indicating the allowable proportion of NA 
#' values within each moving window. If the proportion of NA values exceeds this, the window's value is set as 
#' NA; otherwise, the computation uses the non-NA values.
#' @param rescale Logical; if TRUE, scales and centres the values in each element of 'x'.
#' @param lambda Numeric; lambda value used for Minkowski distance calculation.
#' @param diag Logical; if TRUE, includes the diagonal of the distance matrix in computations.
#' @param time_vector time; 
#' @param steepness numeric; steepness of the logistic function.
#' @param midpoint numeric; midpoint of the logistic function
#' @param cycle_length string; The length of the cycle. Can be a numeric value or a string specifying the units ('year', 'month', 'day', 'hour', 'minute', 'second'). When numeric, the cycle length is in the same units as time_scale. When a string, it specifies the time unit of the cycle.
#' @param time_scale string; Specifies the time scale for the conversion. Must be one of 'year', 'month', 'day', 'hour', 'minute', 'second'. When cycle_length is a string, time_scale changes the unit in which the result is expressed. When cycle_length is numeric, time_scale is used to compute the elapsed time in seconds.
#' @param debugging Logical; if TRUE, additional diagnostic messages are output, useful for debugging. Default 
#' is FALSE.
#' @param isfloat Logical; specifies if the input data are floats.
#' @param mfactor Numeric; multiplication factor applied if input data are float numbers.
#' @param np Number of processes for parallel computation.
#' @param progBar logical. If TRUE a progress bar is shown.
#' @return A list of matrices, each representing a layer of the input RasterStack, containing calculated 
#' Rao's index values. The dimensions correspond to those of the input, and the list length is equal to the 
#' length of 'alpha'.
#' @seealso \code{\link{paRao}} for the parallelized version of the Rao's index computation.
#' @author Duccio Rocchini \email{duccio.rocchini@@unibo.it}, 
#' Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}

mpaRaoS <- function(x, alpha, window, dist_m, na.tolerance, rescale, lambda,
                    diag, time_vector, stepness, midpoint, cycle_length,
                    time_scale, debugging, isfloat, mfactor, np, progBar) {

  win   <- window
  NAwin <- 2 * window + 1
  full_window_n <- NAwin * NAwin

  message("\n\nProcessing alpha: ", alpha, " Moving Window: ", NAwin)

  if (progBar) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent in :elapsed\n",
      total = ncol(x[[1]]),
      clear = FALSE,
      width = 60,
      force = FALSE
    )
  }

  mfactor <- ifelse(isfloat, mfactor, 1)

  # use first layer only for dimensions
  rasterm <- x[[1]]
  nrow_x <- nrow(rasterm)
  ncol_x <- ncol(rasterm)

  # Validate distance function
  validDistanceMetrics <- c("euclidean", "manhattan", "canberra",
                            "minkowski", "mahalanobis", "twdtw")

  if (is.character(dist_m) && dist_m %in% validDistanceMetrics) {
    distancef <- switch(
      dist_m,
      euclidean   = get(".meuclidean"),
      manhattan   = get(".mmanhattan"),
      canberra    = get(".mcanberra"),
      minkowski   = {
        if (lambda == 0) {
          stop("Minkowski distance with lambda = 0 is undefined. Choose another value.")
        }
        get(".mminkowski")
      },
      mahalanobis = {
        warning("Mahalanobis distance is not fully supported for multidimensional Rao's Q.")
        get(".mmahalanobis")
      },
      twdtw       = get(".mtwdtw")
    )
  } else if (is.matrix(dist_m)) {
    distancef <- dist_m
  } else {
    stop("Invalid distance metric. Choose among 'euclidean', 'manhattan', 'canberra', 'minkowski', 'mahalanobis', 'twdtw', or provide a matrix.")
  }

  if (debugging) {
    message("#check: distance function set.")
  }

  # Pad all layers once
  hor <- matrix(NA, ncol = ncol_x, nrow = win)
  ver <- matrix(NA, ncol = win, nrow = nrow_x + win * 2)

  trastersm <- lapply(x, function(layer) {
    cbind(ver, rbind(hor, layer, hor), ver)
  })

  if (debugging) {
    message("#check: padded layers built.")
  }

  # Output
  raoqe <- matrix(NA_real_, nrow = nrow_x, ncol = ncol_x)

  # minimum number of non-NA cells required in a window
  min_valid_cells <- floor(full_window_n - (full_window_n * na.tolerance))

  # helper to aggregate distances without eval(parse())
  aggregate_alpha <- function(vout, alpha, full_window_n) {
    if (length(vout) == 0 || all(is.na(vout))) {
      return(NA_real_)
    }

    if (alpha >= .Machine$integer.max || is.infinite(alpha)) {
      return(max(vout * 2, na.rm = TRUE))
    }

    if (alpha > 0) {
      if (alpha > 100) {
        warning("With this alpha value you may get integer overflow. Consider decreasing the value of alpha.")
      }
      return((sum(rep(vout^alpha, 2) * (1 / (full_window_n^2)), na.rm = TRUE))^(1 / alpha))
    }

    if (alpha == 0) {
      return(prod(vout, na.rm = TRUE)^(1 / (full_window_n^2)))
    }

    stop("alpha must be >= 0.")
  }

  # helper for TWDTW distance
  twdtw_pair <- function(a, b) {
    distancef(
      x = list(a, b),
      time_vector = time_vector,
      stepness = stepness,
      midpoint = midpoint,
      cycle_length = cycle_length,
      time_scale = time_scale
    ) / mfactor
  }

  # helper for non-TWDTW distance
  other_pair <- function(a, b) {
  lpair <- lapply(seq_along(a), function(k) c(a[k], b[k]))
  distancef(lpair) / mfactor
}

  total_cells <- ncol_x * nrow_x * length(x)
  if (total_cells > 10000) {
    message("\n Warning: ", total_cells, " cells to be processed, it may take some time...\n")
  }

  for (cl in (1 + win):(ncol_x + win)) {
    if (progBar) pb$tick()

    for (rw in (1 + win):(nrow_x + win)) {

      # extract local windows from all layers
      tw <- lapply(trastersm, function(layer) {
        layer[(rw - win):(rw + win), (cl - win):(cl + win), drop = FALSE]
      })

      if (dist_m == "twdtw") {
        # For TWDTW, use only complete trajectories across all time steps
        traj_mat <- do.call(
          cbind,
          lapply(tw, function(m) as.vector(t(m)))
        )

        complete_idx <- stats::complete.cases(traj_mat)
        traj_mat <- traj_mat[complete_idx, , drop = FALSE]

        if (nrow(traj_mat) < min_valid_cells) {
          raoqe[rw - win, cl - win] <- NA_real_
          next
        }

        if (nrow(traj_mat) < 2) {
          raoqe[rw - win, cl - win] <- 0
          next
        }

        # Optional short-circuit: all trajectories identical
        if (nrow(unique(traj_mat)) < 2) {
          raoqe[rw - win, cl - win] <- 0
          next
        }

        ntraj <- nrow(traj_mat)
        npairs <- ntraj * (ntraj - 1) / 2
        vout <- numeric(npairs)

        k <- 1L
        for (i in 1:(ntraj - 1L)) {
          ai <- traj_mat[i, ]
          for (j in (i + 1L):ntraj) {
            vout[k] <- twdtw_pair(ai, traj_mat[j, ])
            k <- k + 1L
          }
        }

        raoqe[rw - win, cl - win] <- aggregate_alpha(vout, alpha, full_window_n)

      } else {
        # Keep the non-TWDTW branch close to the old behavior
        # but still use one trajectory matrix per window
        traj_mat <- do.call(
          cbind,
          lapply(tw, function(m) as.vector(t(m)))
        )

        # old code effectively tolerated partial missingness poorly;
        # here we exclude incomplete rows for stability
        complete_idx <- stats::complete.cases(traj_mat)
        traj_mat <- traj_mat[complete_idx, , drop = FALSE]

        if (nrow(traj_mat) < min_valid_cells) {
          raoqe[rw - win, cl - win] <- NA_real_
          next
        }

        if (nrow(traj_mat) < 2) {
          raoqe[rw - win, cl - win] <- 0
          next
        }

        if (nrow(unique(traj_mat)) < 2) {
          raoqe[rw - win, cl - win] <- 0
          next
        }

        ntraj <- nrow(traj_mat)
        npairs <- ntraj * (ntraj - 1) / 2
        vout <- numeric(npairs)

        k <- 1L
        for (i in 1:(ntraj - 1L)) {
          ai <- traj_mat[i, ]
          for (j in (i + 1L):ntraj) {
            vout[k] <- other_pair(ai, traj_mat[j, ])
            k <- k + 1L
          }
        }

        raoqe[rw - win, cl - win] <- aggregate_alpha(vout, alpha, full_window_n)
      }
    }
  }

  return(raoqe)
}