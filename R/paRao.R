
#' Parametric Rao's index of quadratic entropy (Q)
#'
#' It computes the parametric version of Rao's index of quadratic entropy (Q) on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects.
#' @param area Input vector area layer for area-based calculation.
#' @param field Column name of the vector area layer to use to calculate the index.
#' @param dist_m Define the type of distance to be calculated between numerical categories. `dist_m` can be a character string which defines the name of the distance to derive such as "euclidean". The distance names allowed are the same as for \code{proxy::dist}. Alternatively, `dist_m` can be a function which calculates a user-defined distance, (i.e., \code{function(x,y) {return(cos(y-x)-sin(y-x))}}) or a matrix of distances. If `method="multidimension"` then only "euclidean", "manhattan", "canberra", "minkowski" and "mahalanobis" can be used. Default value is "euclidean". If `dist_m` is a matrix, then the function will assume that the matrix contains the distances. Moreover \emph{"twdtw"} (time weighted dynamic time warping) can be used as a way to calculate distances for time series in the `paRao` multidimensional mode. 
#' @param window The side of the square moving window, it must be a vector of odd numeric values greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3. `window` can be a vector with length greater than 1, in this case, Rao's index will be calculated over `x` for each value in the vector.
#' @param alpha Weight for the distance matrix. If `alpha = 0`, distances will be averaged with a geometric average, if `alpha=1` with an arithmetic mean, if `alpha = 2` with a quadratic mean, `alpha = 3` with a cubic mean, and so on. if `alpha` tends to infinite (i.e., higher than the maximum integer allowed in R) or `alpha=Inf`, then the maximum distance will be taken. `alpha` can be a vector with length greater than 1, in this case, Rao's index will be calculated over `x` for each value in the vector.
#' @param method Currently, there are two ways to calculate the parametric version of Rao's index. If `method="classic"`, then the normal parametric Rao's index will be calculated on a single matrix. If `method="multidimension"` (experimental!), a list of matrices must be provided as input. In the latter case, the overall distance matrix will be calculated in a multi- or hyper-dimensional system by using the distance measure defined through the function argument `dist_m`. Each pairwise distance is then multiplied by the inverse of the squared number of pixels in the considered moving window, and the Rao's Q is finally derived by applying a summation. Default value is `"classic"`.
#' @param rasterOut Boolean, if TRUE the output will be a SpatRaster object with `x` as a template.
#' @param lambda The value of the lambda of Minkowski's distance. Considered only if `dist_m = "minkowski"` and `method="multidimension"`. Default value is 0.
#' @param na.tolerance Numeric value (0.0-1.0) which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over `x`. If the relative proportion of NA's in a moving window is bigger than `na.tolerance`, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default value is 1.0. In the rare event that a moving window has NA cells number < `na.tolerance` threshold but has only 1 non-NA value, then its resulting Rao value will always be 0.
#' @param rescale Boolean. Considered only if `method="multidimension"`. If TRUE, each element of `x` is rescaled and centred.
#' @param diag Boolean. If TRUE then the diagonal of the distance matrix is filled with 0's, otherwise with NA's. If `diag=TRUE` and `alpha=0`, the output matrix will inexorably be 0's.
#' @param simplify Number of decimal places to be retained to calculate distances in Rao's index. Default `simplify=0`.
#' @param np The number of processes (cores) which will be spawned. Default value is 2.
#' @param cluster.type The type of cluster which will be created. The options are `"MPI"` (which calls "makeMPIcluster"), `"FORK"`, and `"SOCK"` (which call "makeCluster"). Default type is `"SOCK"`.
#' @param progBar logical. If TRUE a progress bar is shown.
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.
#' @param time_vector time; 
#' @param stepness numeric; steepness of the logistic function.
#' @param midpoint numeric; midpoint of the logistic function
#' @param cycle_length string; The length of the cycle. Can be a numeric value or a string specifying the units ('year', 'month', 'day', 'hour', 'minute', 'second'). When numeric, the cycle length is in the same units as time_scale. When a string, it specifies the time unit of the cycle.
#' @param time_scale string; Specifies the time scale for the conversion. Must be one of 'year', 'month', 'day', 'hour', 'minute', 'second'. When cycle_length is a string, time_scale changes the unit in which the result is expressed. When cycle_length is numeric, time_scale is used to compute the elapsed time in seconds.
#' @return A list of matrices of dimension `dim(x)` with length equal to the length of `alpha`. If `rasterOut=TRUE` and `x` is a SpatRaster, then the output is a list of SpatRaster objects.
#' @details The parametric Rao's Index (Q) is an extension of Rao's Index which considers a generalized mean between distances. The general formula for the parametric Rao's index is Q_a = \deqn{Q = \sum_{i, j} p_i p_j d_{ij}^{\alpha}}. Where `N` is the number of numerical categories, `i` and `j` are pair of numerical categories in the same moving window, and `alpha` is a weight given to distances. In the "multidimension" Rao's index, first the distances among categories are calculated considering more than one feature, and then the overall Rao's Q is derived by using these distances.
#' @references 
#' Rao, C. R. (1982). Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology, 21(1), 24-43. 
#' 
#' @examples
#' \dontrun{
#' # loading data
#' data(volcano)
#' r <- terra::rast(volcano)
#' 
#' # we want to compute Rao's index on this data using a 3x3 window
#' res <- paRao(x = r, window = 3, alpha = 2, method = "classic")
#' terra::plot(res[[1]][[1]])
#' }
#'
#' @export

paRao <- function(x, area = NULL, field = NULL, dist_m = "euclidean",
                  window = 9, alpha = 1, method = "classic",
                  rasterOut = TRUE, lambda = 0, na.tolerance = 1.0,
                  rescale = FALSE, diag = TRUE, simplify = 0,
                  np = 1, cluster.type = "SOCK", progBar = TRUE,
                  debugging = FALSE, time_vector = NULL,
                  stepness = -0.5, midpoint = 35,
                  cycle_length = "year", time_scale = "day") {

  isfloat <- FALSE
  israst  <- FALSE
  mfactor <- 1

  method <- match.arg(method, c("classic", "multidimension"))

  if (method == "multidimension") {
    warning("Multidimension Rao's index is experimental and should be used with caution.")
  }

  if (!is.null(simplify) && simplify > 0) {
    warning(paste0("Simplify=", simplify, ". Rounding data to ", simplify, " decimal places."))
  }

  if (!(methods::is(x, "matrix") || methods::is(x, "SpatRaster") || methods::is(x, "list"))) {
    stop("Invalid input: x must be a matrix, a SpatRaster, or a list.")
  }

  if (!is.numeric(na.tolerance) || length(na.tolerance) != 1L || is.na(na.tolerance) ||
      na.tolerance < 0 || na.tolerance > 1) {
    stop("na.tolerance must be a numeric value in the [0, 1] interval.")
  }

  if (!is.numeric(alpha) || any(is.na(alpha))) {
    stop("alpha must be a numeric vector with no NA values.")
  }
  if (any(alpha < 0)) {
    stop("alpha values must be non-negative.")
  }

  if (!is.null(area)) {
    if (!methods::is(area, "SpatVector")) {
      stop("area must be a SpatVector.")
    }
    if (is.null(field) || !field %in% names(area)) {
      stop("field must be a valid variable name in 'area'.")
    }
    if (np > 1) {
      stop("Parallel area-based Rao's index is not yet implemented.")
    }
    message("Processing area-based Rao's index.")
  }

  # Normalize input structure
  if (method == "classic") {
    if (methods::is(x, "list")) {
      stop("For list input, method must be 'multidimension'.")
    }
    rasterm <- list(x)
  } else {
    if (methods::is(x, "list")) {
      rasterm <- x
    } else if (methods::is(x, "SpatRaster")) {
      rasterm <- lapply(seq_len(terra::nlyr(x)), function(i) x[[i]])
    } else {
      stop("For method = 'multidimension', x must be a list or a multi-layer SpatRaster.")
    }
  }

  # TWDTW checks
  if (dist_m == "twdtw") {
    if (method != "multidimension") {
      stop("dist_m = 'twdtw' requires method = 'multidimension'.")
    }
    if (!is.list(rasterm)) {
      stop("For dist_m = 'twdtw', x must be a list of time-ordered layers.")
    }
    if (is.null(time_vector)) {
      stop("time_vector must be provided if dist_m = 'twdtw'.")
    }
    if (length(time_vector) != length(rasterm)) {
      stop("time_vector must have the same length as the number of layers in x.")
    }
  }

  # Window validation via helper
  w <- calculateWindow(window)

  # Convert input layers to matrices
  if (is.null(area)) {
    if (all(sapply(rasterm, methods::is, "SpatRaster"))) {
      israst <- TRUE

      if (!all(sapply(rasterm, terra::is.int))) {
        warning("Input data are float numbers. Converting data to integer matrices.")
        isfloat <- TRUE
        mfactor <- 100^simplify
      }

      rasterm <- lapply(rasterm, function(z) {
        if (rescale && method == "multidimension") {
          message("Centring and scaling data...")
          z <- terra::scale(z, center = TRUE, scale = TRUE)
        }

        y <- terra::as.matrix(z, wide = TRUE)

        if (!all(y == round(y), na.rm = TRUE)) {
          isfloat <<- TRUE
          mfactor <<- 100^simplify
          y <- round(y * mfactor)
        }

        storage.mode(y) <- "integer"
        y
      })

    } else if (all(sapply(rasterm, methods::is, "matrix"))) {
      israst <- FALSE

      if (!all(sapply(rasterm, function(m) all(m == round(m), na.rm = TRUE)))) {
        warning("Input data are float numbers. Converting data to integer matrices...")
        isfloat <- TRUE
        mfactor <- 100^simplify
      }

      rasterm <- lapply(rasterm, function(z) {
        if (rescale && method == "multidimension") {
          message("Centring and scaling data...")
          z <- (z - mean(z, na.rm = TRUE)) / stats::sd(z, na.rm = TRUE)
        }

        if (!all(z == round(z), na.rm = TRUE)) {
          y <- round(z * mfactor)
        } else {
          y <- z
        }

        storage.mode(y) <- "integer"
        y
      })

    } else {
      stop("All elements of x must be of the same supported class.")
    }
  }

  if (np > 1 && progBar) {
    message("Progress bar disabled for parallel execution.")
  }

  # Run functions
  if (np == 1) {

    if (method == "classic") {
      if (!is.null(area)) {
        if (debugging) cat("#check: Inside classic area clause.")
        split_layers <- terra::split(area, field)
        out <- lapply(split_layers, function(are) {
          lapply(alpha, function(a) {
            paRaoAreaS(
              area = are,
              rasterm = rasterm[[1]],
              simplify = simplify,
              alpha = a
            )
          })
        })
      } else {
        out <- lapply(w, function(win) {
          lapply(alpha, function(a) {
            paRaoS(
              x = rasterm[[1]],
              alpha = a,
              window = win,
              dist_m = dist_m,
              na.tolerance = na.tolerance,
              diag = diag,
              debugging = debugging,
              isfloat = isfloat,
              mfactor = mfactor,
              progBar = progBar
            )
          })
        })
      }

    } else if (method == "multidimension") {
      if (!is.null(area)) {
        if (debugging) cat("#check: Inside multi area clause.")
        split_layers <- terra::split(area, field)
        out <- lapply(split_layers, function(are) {
          lapply(alpha, function(a) {
            mpaRaoAreaS(
              area = are,
              rasterm = rasterm,
              dist_m = dist_m,
              simplify = simplify,
              alpha = a
            )
          })
        })
      } else {
        out <- lapply(w, function(win) {
          lapply(alpha, function(a) {
            mpaRaoS(
              x = rasterm,
              alpha = a,
              window = win,
              dist_m = dist_m,
              na.tolerance = na.tolerance,
              rescale = rescale,
              lambda = lambda,
              diag = diag,
              time_vector = time_vector,
              stepness = stepness,
              midpoint = midpoint,
              cycle_length = cycle_length,
              time_scale = time_scale,
              debugging = debugging,
              isfloat = isfloat,
              mfactor = mfactor,
              np = np,
              progBar = progBar
            )
          })
        })
      }
    }

  } else {

    cls <- openCluster(cluster.type, np, progBar, debugging)
    on.exit(stopCluster(cls), add = TRUE)
    gc()

    if (method == "classic") {
      out <- lapply(w, function(win) {
        lapply(alpha, function(a) {
          paRaoP(
            x = rasterm[[1]],
            alpha = a,
            window = win,
            dist_m = dist_m,
            na.tolerance = na.tolerance,
            diag = diag,
            debugging = debugging,
            isfloat = isfloat,
            mfactor = mfactor,
            np = np,
            progBar = progBar
          )
        })
      })

    } else if (method == "multidimension") {
      out <- lapply(w, function(win) {
        lapply(alpha, function(a) {
          mpaRaoP(
            x = rasterm,
            alpha = a,
            window = win,
            dist_m = dist_m,
            na.tolerance = na.tolerance,
            rescale = rescale,
            lambda = lambda,
            diag = diag,
            time_vector = time_vector,
            stepness = stepness,
            midpoint = midpoint,
            cycle_length = cycle_length,
            time_scale = time_scale,
            debugging = debugging,
            isfloat = isfloat,
            mfactor = mfactor,
            np = np,
            progBar = progBar
            )
          })
        })
    }
  }

  # Format output
  if (!is.null(area)) {
    y <- do.call(rbind.data.frame, lapply(out, function(x) rbind(x)))
    if (nrow(y) > 1) y <- as.data.frame(sapply(y, unlist))
    names(y) <- paste("alpha.", alpha, sep = "")
    terra::values(area) <- cbind.data.frame(area, y)
    return(area)
  }

  if (rasterOut && israst) {
    outR <- lapply(out, function(insm) {
      if (method == "multidimension") {
        y <- lapply(insm, terra::rast, crs = terra::crs(x[[1]]), ext = terra::ext(x[[1]]))
      } else {
        y <- lapply(insm, terra::rast, crs = terra::crs(x), ext = terra::ext(x))
      }
      names(y) <- paste("alpha.", alpha, sep = "")
      y
    })
    names(outR) <- paste("window.", window, sep = "")
    return(outR)
  }

  outM <- lapply(out, function(insm) {
    names(insm) <- paste("alpha.", alpha, sep = "")
    insm
  })
  names(outM) <- paste("window.", window, sep = "")
  return(outM)
}
#' Rao's index
#'
#' An alias for `paRao` with `alpha` fixed at 2.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects.
#' @param ... Other parameters passed to `paRao`.
#' @return A return value description.
#' @examples
#' \dontrun{
#' data(volcano)
#' r <- terra::rast(volcano)
#' res <- Rao(x = r, window = 3)
#' terra::plot(res[[1]][[1]])
#' }
#' @export
Rao <- function(x, ...) {
  paRao(x = x, alpha = 2, ...)
}
