#' Sequential Shannon's diversity index
#'
#' This function computes Shannon's diversity index using a sequential method. It is designed for situations 
#' where parallel computation is not feasible or desired. The function applies a moving window approach to 
#' the provided raster data.
#'
#' @param rasterm Input data; typically a matrix of raster data values.
#' @param w Numeric; half of the side of the square moving window. It determines the area over which each 
#' Shannon index value is calculated.
#' @param na.tolerance Numeric; a tolerance threshold (between 0.0 and 1.0) for NA values in the moving 
#' window. If the proportion of NA values in a window exceeds this threshold, the result for that 
#' window is set as NA. Otherwise, the calculation ignores the NA values. The default is 0.0, 
#' indicating no tolerance for NA values.
#' @param debugging Boolean; if TRUE, the function outputs additional information useful for debugging 
#' purposes. This parameter is FALSE by default, meaning that the extra information is not displayed 
#' during normal operation.
#' @return A matrix or a list of matrices, each containing the Shannon diversity index values
#' calculated over the corresponding area of the input raster data.
#' @seealso \code{\link{Shannon}} for the standard (non-sequential) version of the Shannon diversity index calculation.
#' @author Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}, 
#' Martina Iannacito \email{martina.iannacito@@inria.fr}, 
#' Duccio Rocchini \email{duccio.rocchini@@unibo.it}

ShannonS <- function(rasterm, w, na.tolerance, debugging){
  # Reshape values
  out <- matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  values <- as.numeric(as.factor(rasterm))
  rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional  columns and rows for moving window
  #
  hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm <- cbind(ver,rbind(hor,rasterm_1,hor),ver)
  window = 2*w+1
  #
  ## Loop over all the pixels
  #
  for (cl in (1+w):(ncol(rasterm)+w)) {
    for(rw in (1+w):(nrow(rasterm)+w)) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        out[rw-w,cl-w]<-NA
      } else {
        tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if(debugging) {
          message("Shannon-Wiener\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window)
        }
        tw_values<-as.vector(tw)
        p<-tw_values/sum(tw_values)
        out[rw-w,cl-w]<-(-(sum(p*log(p))))
      }
    }
    svMisc::progress(value=cl/(ncol(trasterm)-1)*100, max.value=100, progress.bar = F,init=T)
  } 
  
  return(out)
}