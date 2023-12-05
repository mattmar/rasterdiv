#' Sequential Hill's diversity index
#'
#' Computes Hill's diversity index in a non-parallelized (sequential) manner.
#'
#' @param x Input data.
#' @param window Half of the side of the square moving window.
#' @param alpha Alpha value for the order of diversity in Hill's Index.
#' @param na.tolerance A numeric value between 0.0 and 1.0, which indicates the proportion of NA values that will be tolerated to calculate Hill's index in each moving window over \code{x}. If the relative proportion of NA's in a moving window is greater than na.tolerance, then the value of the window will be set as NA; otherwise, the index will be calculated considering the non-NA values. Default value is 1.0 (i.e., full tolerance for NA's).
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed for debugging purposes.
#'
#' @return Matrix or a list of matrices with the Hill index computed through a moving window of the given size.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@@gmail.com}, Martina Iannacito \email{martina.iannacito@@inria.fr}, Duccio Rocchini \email{duccio.rocchini@@unibo.it}
#'
#' @seealso \code{\link{Hill}}
#'
#' @keywords internal

HillS <- function(x, window = 1, alpha=1, na.tolerance=1, debugging=FALSE){
  # `win` is the operative moving window
  win = window 
  NAwin <- 2*window+1
  message("\n\nProcessing alpha: ",alpha, " Moving Window: ", NAwin)
  # Set a progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed",
    total = (dim(x)[2]+NAwin), 
    clear = FALSE, 
    width = 60, 
    force = FALSE)

  out <- matrix(rep(NA,dim(x)[1]*dim(x)[2]),nrow=dim(x)[1],ncol=dim(x)[2])
  # Reshape values
  values <- as.numeric(as.factor(x))
  x_1 <- matrix(data=values,nrow=dim(x)[1],ncol=dim(x)[2])
  # Add additional columns and rows for moving window
  hor <- matrix(NA,ncol=dim(x)[2],nrow=win)
  ver <- matrix(NA,ncol=win,nrow=dim(x)[1]+win*2)
  tx <- cbind(ver,rbind(hor,x_1,hor),ver)
  # Loop over each pixel
  for (cl in (1+win):(dim(x)[2]+win)) {
    # Update progress bar
    pb$tick()
    for(rw in (1+win):(dim(x)[1]+win)) {
      if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
        out[rw-win,cl-win]<-NA
      } else {
        tw <- summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if(debugging) {
          message("Hill\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        out[rw-win,cl-win] <- drop(sum(p^alpha))^(1/(1-alpha)) 
      }
    }
  }
  return(out)
}