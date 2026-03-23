#' Multidimensional parallel Parametric Rao's index of quadratic entropy (Q)
#'
#' Multidimensional parametric Rao's index of quadratic entropy (Q).
#'
#' @param x input list.
#' @param alpha alpha value for order of diversity in Hill's Index.
#' @param window half of the side of the square moving window.
#' @param dist_m Type of distance used.
#' @param na.tolerance a numeric value \eqn{(0.0-1.0)} which indicates the proportion
#'   of NA values that will be tolerated to calculate Rao's index in each moving
#'   window over \emph{x}. If the relative proportion of NA's in a moving window is
#'   bigger than na.tolerance, then the value of the window will be set as NA,
#'   otherwise Rao's index will be calculated considering the non-NA values.
#'   Default values is 0.0 (i.e., no tolerance for NA's).
#' @param rescale Scale and centre values in each of the element of x.
#' @param lambda Lambda value for Minkowski distance.
#' @param diag Boolean. Diagonal of the distance matrix.
#' @param time_vector time; 
#' @param stepness numeric; steepness of the logistic function.
#' @param midpoint numeric; midpoint of the logistic function
#' @param cycle_length string; The length of the cycle. Can be a numeric value or a string specifying the units ('year', 'month', 'day', 'hour', 'minute', 'second'). When numeric, the cycle length is in the same units as time_scale. When a string, it specifies the time unit of the cycle.
#' @param time_scale string; Specifies the time scale for the conversion. Must be one of 'year', 'month', 'day', 'hour', 'minute', 'second'. When cycle_length is a string, time_scale changes the unit in which the result is expressed. When cycle_length is numeric, time_scale is used to compute the elapsed time in seconds.
#' @param debugging a boolean variable set to FALSE by default. If TRUE, additional
#'   messages will be printed. For de-bugging only.
#' @param isfloat Are the input data floats?
#' @param mfactor Multiplication factor in case of input data as float numbers.
#' @param np the number of processes (cores) which will be spawned.
#' @param progBar logical. If TRUE a progress bar is shown.
#'
#' @return A list of matrices of dimension \code{dim(x)} with length equal to the
#'   length of \code{alpha}.
#'
#' @author Duccio Rocchini \email{duccio.rocchini@unibo.it}, Marcantonio Matteo
#'   \email{marcantoniomatteo@gmail.com}
#'
#' @seealso \code{\link{paRao}}
#'
#' @keywords internal

mpaRaoP <- function(x, alpha, window, dist_m, na.tolerance, rescale, lambda,
                    diag, time_vector, stepness, midpoint, cycle_length,
                    time_scale, debugging, isfloat, mfactor, np, progBar) {

  win <- window
  NAwin <- 2 * window + 1
  full_window_n <- NAwin * NAwin

  message("\n\nProcessing alpha: ", alpha, " Moving Window: ", NAwin)

  # No worker-side progress bar updates
  if (np > 1 && progBar) {
    progBar <- FALSE
  }

  mfactor <- ifelse(isfloat, mfactor, 1)

  rasterm <- x[[1]]
  nrow_x <- nrow(rasterm)
  ncol_x <- ncol(rasterm)

  validDistanceMetrics <- c(
    "euclidean", "manhattan", "canberra",
    "minkowski", "mahalanobis", "twdtw"
  )

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

  # Pad layers once
  hor <- matrix(NA, ncol = ncol_x, nrow = win)
  ver <- matrix(NA, ncol = win, nrow = nrow_x + win * 2)

  trastersm <- lapply(x, function(layer) {
    cbind(ver, rbind(hor, layer, hor), ver)
  })

  if (debugging) {
    message("#check: padded layers built.")
  }

  # minimum number of valid trajectories required
  min_valid_cells <- floor(full_window_n - (full_window_n * na.tolerance))

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

  other_pair <- function(a, b) {
    lpair <- lapply(seq_along(a), function(k) c(a[k], b[k]))
    distancef(lpair) / mfactor
  }

  total_cells <- ncol_x * nrow_x * length(x)
  if (total_cells > 10000) {
    message("\nWarning: ", total_cells, " cells to be processed, it may take some time...\n")
  }

  out <- foreach::foreach(
    cl = (1 + win):(ncol_x + win),
    .verbose = FALSE
  ) %dopar% {

    col_out <- sapply((1 + win):(nrow_x + win), function(rw) {

      # extract local windows from all layers
      tw <- lapply(trastersm, function(layer) {
        layer[(rw - win):(rw + win), (cl - win):(cl + win), drop = FALSE]
      })

      # build trajectory matrix: rows = local pixels, cols = layers/time
      traj_mat <- do.call(
        cbind,
        lapply(tw, function(m) as.vector(t(m)))
      )

      # complete trajectories only
      traj_mat <- traj_mat[stats::complete.cases(traj_mat), , drop = FALSE]

      # NA tolerance check
      if (nrow(traj_mat) < min_valid_cells) {
        return(NA_real_)
      }

      # fewer than 2 valid trajectories
      if (nrow(traj_mat) < 2) {
        return(0)
      }

      # all trajectories identical
      if (nrow(unique(traj_mat)) < 2) {
        return(0)
      }

      ntraj <- nrow(traj_mat)
      npairs <- ntraj * (ntraj - 1) / 2
      vout <- numeric(npairs)

      k <- 1L
      for (i in 1:(ntraj - 1L)) {
        ai <- traj_mat[i, ]

        for (j in (i + 1L):ntraj) {
          bj <- traj_mat[j, ]

          if (dist_m == "twdtw") {
            vout[k] <- twdtw_pair(ai, bj)
          } else {
            vout[k] <- other_pair(ai, bj)
          }

          k <- k + 1L
        }
      }

      aggregate_alpha(vout, alpha, full_window_n)
    })

    col_out
  }

  do.call(cbind, out)
}