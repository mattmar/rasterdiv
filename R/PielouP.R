#' Parallelised Pielou's diversity index
#'
#' This function calculates Pielou's diversity index in a parallelized manner, 
#' allowing for improved performance on suitable hardware. The diversity index 
#' is computed using a moving window approach over the input data.
#'
#' @param rasterm Input raster data, representing the environmental variable(s) 
#' over which the diversity index should be calculated.
#' @param w The size of the half-side of the square moving window used in the 
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
#' @seealso \code{\link{Pielou}} for the non-parallelized version of the 
#' Pielou's diversity index calculation.
#'
#' @examples
#' \dontrun{
#' # Demonstration of function with hypothetical data
#' # Ensure you replace this with actual raster data
#' demo_raster <- #... (your raster data here)
#' result <- PielouP(rasterm = demo_raster, w = 3, na.tolerance = 0.1, debugging = FALSE)
#' # proceed with analyzing 'result'
#' }

PielouP<-function(rasterm, w, na.tolerance, debugging){
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
  
  PielouOP <- foreach::foreach(cl=(1+w):(ncol(rasterm)+w),.options.parallel = opts,.verbose = FALSE) %dopar% {
    PielouOut <- sapply((1+w):(nrow(rasterm)+w), function(rw) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        vv<-NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }

        if( debugging ) {
          message("\nPielou - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
        }

        tw_values <- as.vector(tw)
        maxS <- log(length(tw))
        p <- tw_values/sum(tw_values)
        vv <- (-(sum(p*log(p))))/maxS

        if( debugging ) {
          message("\ncat: ",paste(names(tw),collapse=" ")," log S: ",maxS," Pielou: ",vv)
        }

        return(vv)
      }
    })
    return(PielouOut)
  }
  return(PielouOP)
}