#' Parallelised Renyi's diversity index
#'
#' This function provides a parallelised computation of Renyi's diversity index.
#'
#' @param rasterm Input data, expected to be in a specific format suitable for the calculations.
#' @param w Half of the side of the square moving window used in the computation.
#' @param alpha The alpha value for the order of diversity in Renyi's Index calculations.
#' @param base The base of the logarithm used in the diversity index calculations.
#' @param na.tolerance A numeric value between 0.0 and 1.0 indicating the proportion of NA values 
#' that are acceptable in each moving window over \code{rasterm} during the calculation of Rao's index.
#' If the proportion of NAs in a window exceeds this value, the result for that window is set as NA.
#' Otherwise, the index is calculated using the non-NA values. The default is 0.0, indicating no tolerance for NAs.
#' @param debugging Boolean; controls the generation of additional, diagnostic messages. 
#' If TRUE, extra messages are printed primarily for debugging purposes. Default is FALSE.
#' @return A matrix or a list of matrices, each containing the Renyi index values computed using a moving 
#' window of the specified size.
#' @author Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}, 
#' Martina Iannacito \email{martina.iannacito@@inria.fr}, 
#' Duccio Rocchini \email{duccio.rocchini@@unibo.it}
#' @seealso \code{\link{Renyi}} for the non-parallelised version of the Renyi diversity index calculation.
#' @keywords internal

RenyiP <- function(rasterm, w, alpha, base, na.tolerance, debugging){
  # Some initial housekeeping
  window = 2*w+1
  message("\n\nProcessing alpha ",alpha, " Window ", window)
  # Set a progress bar
  pb <- utils::txtProgressBar(title = "Iterative training", min = w, max = dim(rasterm)[2]+w, style = 3)
  #
  ## Reshape values
  #
  values <- as.numeric(as.factor(rasterm))
  rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm <- cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Start the parallelized loop over iter
  #
  RenyiOP <- foreach::foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    utils::setTxtProgressBar(pb, cl)
    if(debugging) {
      cat(paste(cl))
    }
    RenyiOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        vv <- NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if( debugging ) {
          message("Renyi - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        vv <- 1/(1-alpha) * drop(log(sum(p^alpha),base))
        return(vv)
      }
    })
    return(RenyiOut)
  } # End Renyi - parallelised
  message(("\n\n Parallel calculation of Renyi's index complete.\n"))
  return(RenyiOP)
}