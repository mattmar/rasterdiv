#' Sequential Berger-Parker's diversity index
#'
#' Berger-Parker's diversity index calculated sequentially over a raster matrix.
#'
#' @param rasterm Input data in the form of a raster matrix.
#' @param w The half of the side of the square moving window. It defines the size of the area (neighbourhood) over which the diversity index is calculated.
#' @param na.tolerance A numeric value between 0.0 and 1.0, indicating the proportion of NA values that will be tolerated in each moving window over \emph{rasterm} during the calculation of the diversity index. If the proportion of NA's in a moving window exceeds na.tolerance, the value for that window will be set as NA; otherwise, the diversity index is calculated considering the non-NA values. The default value is 0.0, indicating no tolerance for NA's.
#' @param debugging A boolean that controls the display of additional messages for debugging purposes. Set to FALSE by default. If TRUE, messages detailing the function's progress and computations will be printed.
#'
#' @return A matrix or a list of matrices containing the Berger-Parker diversity index values, computed for each pixel using a moving window of the specified size.
#'
#' @seealso \code{\link{BergerParker}} for the non-sequential version of the Berger-Parker diversity index calculation.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@gmail.com}, Martina Iannacito \email{martina.iannacito@inria.fr}, Duccio Rocchini \email{duccio.rocchini@unibo.it}
#' 
#' @examples
#' \dontrun{
#' # Create a sample raster matrix
#' a <- matrix(c(10, 10, 10, 20, 20, 20, 20, 30, 30), ncol = 3, nrow = 3)
#' 
#' # Calculate the Berger-Parker index using a moving window of size 1 (w = 1)
#' result <- BergerParkerS(rasterm = a, w = 1, na.tolerance = 0.1, debugging = FALSE)
#' print(result)
#' }
#' 
#' @importFrom methods is
#' @importFrom stats setNames

BergerParkerS <- function(rasterm, w, na.tolerance, debugging){

  # Reshape values
  out<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  window = 2*w+1
  #
  ## Loop over all the pixels
  #
  for (cl in (1+w):(dim(rasterm)[2]+w)) {
    for(rw in (1+w):(dim(rasterm)[1]+w)) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        out[rw-w,cl-w]<-NA
      } else {
        tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]))
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if(debugging) {
          message("Berger-Parker\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*w+1)
        }
        tw_values<-as.vector(tw)
        out[rw-w,cl-w]<-max(tw_values/sum(tw_values))
      }
    }
    svMisc::progress(value=cl/(ncol(trasterm)-1)*100, max.value=100, progress.bar = F,init=T)
  } 
  return(out)
}