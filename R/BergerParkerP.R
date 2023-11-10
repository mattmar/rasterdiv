#' Parallelised Berger-Parker's diversity index
#'
#' @name BergerParkerP
#' @title Parallelised Berger-Parker's diversity index
#' @description This function computes the Berger-Parker diversity index in a parallelised manner, 
#' improving performance on suitable hardware for large datasets.
#'
#' @param rasterm Input raster data as a matrix or suitable object convertible to a matrix, containing the values on which to compute the index.
#' @param w Integer, representing half of the side length of the square moving window used in calculations.
#' @param na.tolerance Numeric value between 0.0 and 1.0 indicating the proportion of NA values tolerated in each moving window over the input rasterm. If the proportion of NA's exceeds na.tolerance, the result for that window is set as NA. The default is 0.0, indicating no tolerance for NA values.
#' @param debugging Logical, when TRUE, additional console messages are printed for debugging purposes. Defaults to FALSE.
#'
#' @return A matrix or a list of matrices with the Berger-Parker index computed through a moving window of the specified size.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@gmail.com}
#' @author Martina Iannacito \email{martina.iannacito@inria.fr}
#' @author Duccio Rocchini \email{duccio.rocchini@unibo.it}
#'
#' @seealso \code{\link{BergerParker}} for the non-parallel version of the Berger-Parker index calculation.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Generate a sample matrix
#' data_matrix <- matrix(sample(c(NA, 1:10), 100, replace = TRUE), 10, 10)
#'
#' # Compute the Berger-Parker index with parallel processing
#' result <- BergerParkerP(rasterm = data_matrix, w = 1, na.tolerance = 0.1, debugging = TRUE)
#' }
#'

BergerParkerP<-function(rasterm, w, na.tolerance, debugging){
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
  ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  window = 2*w+1
  rm(hor,ver,rasterm_1,values); gc()
  #
  ## Progression bar
  #
  pb <- utils::txtProgressBar(min = (1+w), max = dim(rasterm)[2], style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #
  ## Start the parallelized loop over iter
  #
  BergerParkerOP <- foreach::foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.parallel = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    BergerParkerOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
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
          message("Berger-Parker - parallelised\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        vv <- max(tw_values/sum(tw_values))
      }
    })
    return(BergerParkerOut)
  } # End Berger-Parker - parallelized
  return(BergerParkerOP)
}