#' Parallelised Shannon's diversity index
#'
#' This function computes Shannon's diversity index in a parallelized manner, allowing for
#' enhanced performance on suitable hardware. It is especially useful for large datasets.
#'
#' @param rasterm Input data; a matrix of raster data values.
#' @param w Numeric; half of the side of the square moving window used in the calculation.
#' @param na.tolerance Numeric; a tolerance level (between 0.0 and 1.0) for NA values in the moving
#' window. If the proportion of NA values in a window exceeds this level, the result for that
#' window is set as NA; otherwise, the calculation excludes the NA values. The default is 0.0,
#' which allows no NA values.
#' @param debugging Boolean; if TRUE, additional debugging information is printed during the 
#' function's execution. This is helpful for debugging and is FALSE by default.
#' @return A matrix or a list of matrices, each containing the Shannon diversity index values
#' calculated using a moving window approach.
#' @seealso \code{\link{Shannon}} for the non-parallelized version of the Shannon diversity index.
#' @author Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}, 
#' Martina Iannacito \email{martina.iannacito@@inria.fr}, 
#' Duccio Rocchini \email{duccio.rocchini@@unibo.it}

ShannonP <- function(rasterm, w, na.tolerance, debugging){
  #
  ## Reshape values
  #
  values <- as.numeric( as.factor(rasterm) )
  rasterm_1 <- matrix(data = values, nrow = nrow(rasterm), ncol = ncol(rasterm))
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA, ncol = ncol(rasterm), nrow = w)
  ver <- matrix(NA, ncol = w, nrow = nrow(rasterm)+ w * 2)
  trasterm <- cbind(ver, rbind(hor,rasterm_1,hor), ver)
  window = 2*w+1
  rm(hor, ver, rasterm_1, values); gc()
  #
  ## Progression bar
  #
  pb <- utils::txtProgressBar(min = (1+w), max = ncol(rasterm), style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  ShannonOP <- foreach::foreach(cl=(1+w):(ncol(rasterm)+w),.options.parallel = opts,.verbose = F) %dopar% {
    ShannonOut <- sapply((1+w):(nrow(rasterm)+w), function(rw) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        vv <- NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw <- tw[-length(tw)]
        }
        if( debugging ) {
          message("Shannon - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        vv <- (-(sum(p*log(p))))
        return(vv)
      }
    })
    return(ShannonOut)
  }
  return(ShannonOP)
}