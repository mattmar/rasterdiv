#' Sequential Hill's diversity index
#'
#' Computes Hill's diversity index in a non-parallelized (sequential) manner.
#'
#' @param rasterm Input data.
#' @param w Half of the side of the square moving window.
#' @param alpha Alpha value for the order of diversity in Hill's Index.
#' @param na.tolerance A numeric value between 0.0 and 1.0, which indicates the proportion of NA values that will be tolerated to calculate Hill's index in each moving window over \code{rasterm}. If the relative proportion of NA's in a moving window is greater than na.tolerance, then the value of the window will be set as NA; otherwise, the index will be calculated considering the non-NA values. Default value is 1.0 (i.e., full tolerance for NA's).
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed for debugging purposes.
#'
#' @return Matrix or a list of matrices with the Hill index computed through a moving window of the given size.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@@gmail.com}, Martina Iannacito \email{martina.iannacito@@inria.fr}, Duccio Rocchini \email{duccio.rocchini@@unibo.it}
#'
#' @seealso \code{\link{Hill}}
#'
#' @keywords internal

HillS <- function(rasterm, w, alpha, na.tolerance, debugging){

  out <- matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Reshape values
  values <- as.numeric(as.factor(rasterm))
  rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  # Add additional columns and rows for moving window
  hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm <- cbind(ver,rbind(hor,rasterm_1,hor),ver)
  window = 2*w+1
  # Loop over each pixel
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        out[rw-w,cl-w]<-NA
      } else {
        tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if(debugging) {
          message("Hill\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        out[rw-w,cl-w] <- drop(sum(p^alpha))^(1/(1-alpha)) 
      }
    }
    svMisc::progress(value=cl/(ncol(trasterm)-1)*100, max.value=100, progress.bar = F,init=T)
  }
  return(out)
}