#' Sequential Renyi's diversity index
#'
#' This function calculates Renyi's diversity index sequentially.
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
#' @seealso \code{\link{Renyi}} for the related function in the context of diversity index calculations.
#' @keywords internal

RenyiS <- function(rasterm, w, alpha, base, na.tolerance, debugging){
  # Some initial housekeeping
  window = 2*w+1
  message("\n\nProcessing alpha ",alpha, " Window ", window)
  # Set a progress bar
  pb <- progress::progress_bar$new(
    format = "\n [:bar] :elapsed -- Approximate ETA: :eta \n",
    total = (dim(rasterm)[2]+w), 
    clear = FALSE, 
    width = 80, 
    force = TRUE)
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Reshape values
  values <- as.numeric(as.factor(rasterm))
  rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Add additional columns and rows for moving window
  hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  # Loop over each pixel
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    # Update progress bar
    pb$tick()
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
     if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
      out[rw-w,cl-w]<-NA
    } else {
      tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
      if( "NA's"%in%names(tw) ) {
        tw <- tw[-length(tw)]
      }
      if(debugging) {
        message("Renyi\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
      }
      tw_labels <- names(tw)
      tw_values <- as.vector(tw)
      p <- tw_values/sum(tw_values)
      out[rw-w,cl-w] <- 1/(1-alpha) * drop(log(sum(p^alpha),base))
    }
  }
}
return(out)
}