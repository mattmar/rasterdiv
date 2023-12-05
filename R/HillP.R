#' Parallelised Hill's diversity index
#'
#' Parallelised computation of Hill's diversity index.
#'
#' @param rasterm Input data.
#' @param window Half of the side of the square moving window.
#' @param alpha Alpha value for the order of diversity in Hill's Index.
#' @param na.tolerance A numeric value between 0.0 and 1.0, which indicates the proportion of NA values that will be tolerated to calculate the index in each moving window over \code{rasterm}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA; otherwise, the index will be calculated considering the non-NA values. Default value is 0.0 (i.e., no tolerance for NA's).
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed for debugging purposes.
#' @param np Number of processes for parallel computation.
#'
#' @return Matrix or a list of matrices with the Hill index computed through a moving window of the given size.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@@gmail.com}, Martina Iannacito \email{martina.iannacito@@inria.fr}, Duccio Rocchini \email{duccio.rocchini@@unibo.it}
#'
#' @seealso \code{\link{Hill}}
#'
#' @keywords internal

HillP<-function(rasterm, window, alpha, na.tolerance,debugging,np){
  # `win` is the operative moving window
  win = window 
  NAwin <- 2*window+1
  message("\n\nProcessing alpha: ",alpha, " Moving Window: ", NAwin)
  # Set a progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed\n",
    # Total number of ticks is the number of column +NA columns divided the number of processor.
    total = (dim(rasterm)[2]/np)+5, 
    clear = FALSE, 
    width = 60, 
    force = FALSE)
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(rasterm))
  rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
  #
  ## Add additional columns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=window)
  ver<-matrix(NA,ncol=window,nrow=dim(rasterm)[1]+window*2)
  trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
  rm(hor,ver,rasterm_1,values); gc()  
  #
  ## Start the parallelized loop over iter
  #
  HillOP <- foreach::foreach(cl=(1+window):(dim(rasterm)[2]+window),.verbose = F) %dopar% {
    # Update progress bar
    pb$tick()
    if(debugging) {
      cat(paste(cl))
    }
    HillOut <- sapply((1+window):(dim(rasterm)[1]+window), function(rw) {
      if( length(!which(!trasterm[c(rw-window):c(rw+window),c(cl-window):c(cl+window)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
        vv<-NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(trasterm[c(rw-window):c(rw+window),c(cl-window):c(cl+window)]),maxsum=10000)
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
  message("\n\n Parallel calculation of Hill's index complete.\n")
  return(matrix(unlist(HillOP), ncol = ncol(rasterm), nrow = nrow(rasterm), byrow=FALSE))
}