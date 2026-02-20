#' Validate Input Parameters for Diversity Index Calculation
#'
#' @description Validates the input parameters for diversity index calculation functions. Checks for valid raster types, window sizes, alpha values, and NA tolerance levels.
#'
#' @param x Raster object to be validated.
#' @param window Size of the moving window for calculations.
#' @param alpha Diversity index parameter, default is 1.
#' @param na.tolerance Proportion of acceptable NA values within the window (range: 0 to 1).
#'
#' @return None. Throws an error if any input is invalid.
#' @export
validateInputs <- function(x, window, alpha=1, na.tolerance) {
  if (!isRasterValid(x) || any(!is.wholenumber(window)) || any(window %% 2 == 0) ||
    any(!is.numeric(alpha)) || any(alpha < 0) || 
    na.tolerance < 0 || na.tolerance > 1) {
    stop("Invalid inputs.")
  }
}

#' Check if the Raster Object is Valid
#'
#' @description Checks if the provided object is a valid raster object for diversity index calculation.
#'
#' @param x Object to check.
#'
#' @return TRUE if valid, FALSE otherwise.
#' @noRd
isRasterValid <- function(x) {
  methods::is(x, "matrix") || methods::is(x, "SpatRaster") || 
  methods::is(x, "list") || methods::is(x, "RasterStack")
}

#' Prepare Raster Object for Processing
#'
#' @description Prepares the raster object for processing based on its type.
#'
#' @param x Raster object to prepare.
#'
#' @return A prepared raster object.
#' @noRd
prepareRaster <- function(x) {
  if (methods::is(x, "matrix")) list(x)
  else if (methods::is(x, "SpatRaster")) list(matrix(terra::values(x), ncol = ncol(x), nrow = nrow(x), byrow=TRUE))
  else if (methods::is(x, "list")) list(x[[1]])
  else stop("Invalid raster input object.")
}
#' Calculate the Operative Moving Window Half-Size (win)
#'
#' @description
#' Converts a user-specified odd window size (e.g., 3, 5, 7) to the operative
#' half-size used by the moving-window core (win = (window-1)/2).
#'
#' @param window Integer. The (odd) window size in cells (e.g., 3, 5, 7).
#' @return Integer. The half-size (win).
#' @noRd
calculateWindow <- function(window) {
  if (length(window) != 1L) {
    stop("`window` must be a single integer (odd).")
  }
  if (!is.numeric(window) || is.na(window)) {
    stop("`window` must be numeric.")
  }
  if (window %% 2 == 0) {
    stop("`window` must be an odd integer (e.g., 3, 5, 7).")
  }
  as.integer((window - 1) / 2)
}


#' Format Output of Diversity Index Calculation
#'
#' @description Formats the output of diversity index calculation functions based on user preferences.
#'
#' @param out Raw output from diversity index calculation functions.
#' @param rasterOut Logical indicating whether output should be in raster format.
#' @param x Original raster object.
#' @param alpha Alpha parameter used in the calculation.
#' @param window Window size used in the calculation (scalar) OR a vector/list of labels (advanced use).
#'
#' @return Formatted output (either a single object or a named list).
#' @noRd
formatOutput <- function(out, rasterOut, x, alpha, window) {

  is_spat <- isTRUE(rasterOut) && methods::is(x, "SpatRaster")

  to_rast <- function(m) {
    terra::rast(m, crs = terra::crs(x), ext = terra::ext(x))
  }

  safe_names <- function(n, prefix, labels = NULL) {
    if (is.null(labels)) labels <- seq_len(n)
    if (is.list(labels)) labels <- unlist(labels, use.names = FALSE)
    if (length(labels) == 0L) labels <- seq_len(n)
    rep_len(paste0(prefix, labels), n)
  }

  # ---- Case 1: out is list-of-lists (e.g. windows x alpha) ----
  if (is.list(out) && length(out) > 0L && is.list(out[[1]])) {

    outFormatted <- lapply(out, function(windowList) {
      outAlpha <- lapply(windowList, function(alphaObj) {
        if (is_spat) to_rast(alphaObj) else alphaObj
      })
      names(outAlpha) <- safe_names(length(outAlpha), "alpha.", alpha)
      outAlpha
    })

    names(outFormatted) <- safe_names(length(outFormatted), "window.", window)

    return(if (length(outFormatted) > 1L) outFormatted else outFormatted[[1]])
  }

  # ---- Case 2: out is list (either varying by window OR by alpha) ----
  outFormatted <- lapply(out, function(obj) {
    if (is_spat) to_rast(obj) else obj
  })

  # Name multi-output lists safely
  if (length(outFormatted) > 1L) {
    if (length(outFormatted) == length(alpha)) {
      names(outFormatted) <- safe_names(length(outFormatted), "alpha.", alpha)
    } else if (length(outFormatted) == length(window)) {
      names(outFormatted) <- safe_names(length(outFormatted), "window.", window)
    } else {
      # fallback: name by sequence, never error on naming
      names(outFormatted) <- safe_names(length(outFormatted), "out.", seq_along(outFormatted))
    }
    return(outFormatted)
  }

  # Single output: attach index name
  idx <- gsub("\\((.*$)", "", deparse(sys.calls()[[sys.nframe() - 1]]))
  if (methods::is(outFormatted[[1]], "matrix")) {
    attr(outFormatted[[1]], "index") <- idx
  } else if (methods::is(outFormatted[[1]], "SpatRaster")) {
    names(outFormatted[[1]]) <- rep_len(idx, terra::nlyr(outFormatted[[1]]))
  }

  outFormatted[[1]]
}
#' Check if a Number is a Whole Number
#'
#' @description Checks if a given number is a whole number within a specified tolerance.
#'
#' @param x Number to check.
#' @param tol Tolerance for the check.
#'
#' @return TRUE if x is a whole number within the specified tolerance, FALSE otherwise.
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

#' Open a Parallel Cluster
#'
#' Opens a parallel cluster for computation, registers it for parallel operations,
#' and ensures its closure on script exit.
#'
#' @param cluster.type A character string specifying the type of cluster.
#'   Accepted values are "SOCK", "FORK", or "MPI".
#' @param np An integer specifying the number of processes to be used in the parallel cluster.
#' @param progBar logical. If TRUE a progress bar is shown.
#' @param debugging logical. For developer use.
#'
#' @return An object representing the parallel cluster.
#'
#' @examples
#' \dontrun{
#'   # Open a SOCK cluster with 4 cores
#'   cls <- openCluster("SOCK", 4)
#'   # Your parallel computation code here
#'   # The cluster will automatically close when the script exits
#' }
#'
#' @export
#'
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
openCluster <- function(cluster.type="SOCK", np=2, progBar=TRUE, debugging=FALSE) {

  startParallelCalculationMessage()
  
  if( debugging ){ cat("#check: Before parallel function.") }     
  if(cluster.type %in% c("SOCK", "FORK")) {
    invisible(utils::capture.output( cls <- invisible(parallel::makeCluster(np,type=cluster.type, outfile=ifelse(progBar, ' ', tempfile())))))
  } 
  else if(cluster.type == "MPI") {
    invisible(utils::capture.output( cls <- parallel::makeCluster(np, outfile=ifelse(progBar, ' ', tempfile()), useXDR = FALSE, methods = FALSE)))
  } 
  else {
    stop("Invalid cluster type specified.")
  }

  # Register the cluster for parallel operations
  doParallel::registerDoParallel(cls)

  return(cls)
}


startParallelCalculationMessage <- function() {
  cat("\n")
  cat("---- Starting Parallel Calculation ----\n")
  cat("  [o]  [o]  [o]\n")
  cat("  Code is now running on multiple cores...\n")
  cat("  Please wait for the process to complete.\n\n")
}

# startParallelCalculationMessage <- function() {
#   cat("\n")
#   cat("▬▬ι═══════ﺤ Starting Parallel Calculation ﺤ═══════ι▬▬\n")
#   cat("  • (•_•) •  • ( •_•)>⌐■-■  • (⌐■_■) •\n")
#   cat("  • Code is now running on multiple cores... •\n")
#   cat("  • Please wait for the process to complete. •\n\n")
# }

#
##Supporting function for CRE
#

#' Cumulative Residual Entropy
#'
#' This function calculates the Cumulative Residual Entropy (CRE) for a given set of values.
#' 
#' @param B A numeric vector or matrix representing the values for which CRE is to be calculated.
#' @param base The base of the logarithm used in the calculation. The default is the natural logarithm (e).
#' @return A numeric value representing the CRE.
#' @examples
#' B <- c(1, 2, 3, 4)
#' .CRE_(B)
#' @export

.CRE_<-function(B,base=exp(1))
{
  #Cumulative Residual Entropy
  P=.Prob(B)
  Pcre=.CumRes(P)
  -sum(Pcre*log(Pcre,base)*.Deltas(P), na.rm=TRUE)
}

#' Calculate Point Probability
#'
#' This function computes the probability of each point in a given vector or matrix.
#' 
#' @param C A numeric vector or matrix.
#' @return A vector of probabilities corresponding to each point in `C`.
#' @examples
#' C <- c(1, 1, 2, 2, 3)
#' .Prob(C)
#' @export

.Prob<-function(C)
{
  #Point probability
  if( is.null(dim(C))){L=length(C)
    } else {L=dim(C)[1]}
    table(C)/L
  }

#' Calculate Differences Among Values
#'
#' This function computes the differences among values of a table, used in probability calculations.
#'
#' @param P A numeric vector or matrix representing probabilities.
#' @param first The starting value for difference calculation.
#' @return A vector or matrix of differences.
#' @examples
#' P <- c(0.2, 0.3, 0.5)
#' .Deltas(P)
#' @export

.Deltas<-function(P, first=0)
{
  #Difference among values of a table.
  #For multidimensional table the product is given 
  if ((length(dim(P)))==1){
    delta=c(first,diff(as.numeric(names(P))))} else {
      delta=1
      for (dim in dimnames(P)){
        delta=outer(delta,c(first,diff(as.numeric(dim))))
      }
    }
    drop(delta)
  }

#' Calculate Cumulative Residual Probability
#'
#' This function computes the cumulative residual probability for a given set of probabilities.
#'
#' @param a A numeric vector or matrix representing probabilities.
#' @return A numeric vector or matrix of cumulative residual probabilities.
#' @examples
#' a <- data.frame(V1= c(0.2, 0.3, 0.5), V2 =c(0.2, 0.3, 0.5))
#' .CumRes(a)
#' @export

.CumRes<-function(a)
{
  #Calculate Cumulative Residual Probability
  D=dim(a)
  if (length(dim(a))==1){
    return( rev( cumsum(rev(a)) )) }
    atemp=a
    for(i in 1:length(D))
    {
      aa=.Rev(.Cumsum(.Rev(atemp,i),i),i)
      atemp=aa
    }
    atemp
  }

#' Additional supporting functions like `.Reorder`, `.Cumsum`, `.Rev`
#'
#' These functions provide utility operations like reordering dimensions, 
#' computing cumulative sums, and reversing order along a specific dimension.
#'
#' @param a, ax Additional parameters specific to each function.
#' @param ax Additional parameters specific to each function.
#' @return Output varies depending on the function.
#' @export

.Reorder<-function(a,ax)
{
  #Reordering dimension required after the use of apply over more than one dimension  
  D=dim(a)
  maxdim=length(D)
  reorder=2:maxdim
  end=reorder[(ax):(maxdim-1)]
  if (maxdim==ax){end=reorder[0]}
  reorder=c(reorder[0:(ax-1)],1,end)
  aperm(a,reorder)
  
}
.Cumsum<-function(a,ax=1)
{
  #Cumulative sum along a specific dimension
  D=dim(a)
  dimen=1:length(D)
  .Reorder(apply(a, dimen[-ax],cumsum),ax)
}

.Rev<-function(a,ax)
{
  #Reverse order along a specific dimension
  D=dim(a)
  dimen=1:length(D)
  .Reorder(apply(a, dimen[-ax], rev),ax)
}   