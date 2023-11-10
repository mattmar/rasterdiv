#' Parallelised Hill's diversity index
#'
#' Parallelised computation of Hill's diversity index.
#'
#' @param rasterm Input data.
#' @param w Half of the side of the square moving window.
#' @param alpha Alpha value for the order of diversity in Hill's Index.
#' @param na.tolerance A numeric value between 0.0 and 1.0, which indicates the proportion of NA values that will be tolerated to calculate the index in each moving window over \code{rasterm}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA; otherwise, the index will be calculated considering the non-NA values. Default value is 0.0 (i.e., no tolerance for NA's).
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed for debugging purposes.
#'
#' @return Matrix or a list of matrices with the Hill index computed through a moving window of the given size.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@@gmail.com}, Martina Iannacito \email{martina.iannacito@@inria.fr}, Duccio Rocchini \email{duccio.rocchini@@unibo.it}
#'
#' @seealso \code{\link{Hill}}
#'
#' @keywords internal

HillP<-function(rasterm, w, alpha, na.tolerance,debugging){
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
  HillOP <- foreach::foreach(cl=(1+w):(dim(rasterm)[2]+w),.options.parallel = opts,.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    HillOut <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
      if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
        vv<-NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw <- tw[-length(tw)]
        }
        if( debugging ) {
          message("Hill - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        vv <- drop(sum(p^alpha))^(1/(1-alpha))
        return(vv)
      }
    })
    return(HillOut)
  }
  return(HillOP)
}