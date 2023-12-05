#' Sequential Berger-Parker's diversity index
#'
#' Berger-Parker's diversity index calculated sequentially over a raster matrix.
#'
#' @description
#' This function calculates the Berger-Parker's diversity index for each cell in a matrix,
#' considering a specified moving window around each cell.
#'
#' @param x A numeric matrix representing the data on which the index is to be calculated.
#' @param window The width of the moving window to consider for each cell. The actual window size 
#'        will be `(2 * window + 1) x (2 * window + 1)`. Default is 1.
#' @param na.tolerance The tolerance level for missing data within the moving window. 
#'        A window will be processed only if the proportion of non-missing data is above this threshold. 
#'        Value should be between 0 and 1. Default is 1.
#' @param debugging Boolean flag to enable or disable debugging messages. Default is FALSE.
#'
#' @return A matrix of the same dimensions as `x`, where each cell contains the 
#'         Berger-Parker's diversity index calculated for the window around the cell.
#'
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' bp_index <- BergerParkerS(data, window = 1)
#'
#' @export

BergerParkerS <- function(x, window = 1, na.tolerance=1, debugging=FALSE){
  # `win` is the operative moving window
  win = window 
  NAwin <- 2*window+1
  message("\n\nProcessing moving Window: ", NAwin)
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed",
    total = (dim(x)[2]+NAwin), 
    clear = FALSE, 
    width = 60, 
    force = FALSE)
  out<-matrix(rep(NA,dim(x)[1]*dim(x)[2]),nrow=dim(x)[1],ncol=dim(x)[2])
  values<-as.numeric(as.factor(x))
  x_1<-matrix(data=values,nrow=dim(x)[1],ncol=dim(x)[2])
  #
  ## Add additional columns and rows for moving window
  #
  hor<-matrix(NA,ncol=dim(x)[2],nrow=win)
  ver<-matrix(NA,ncol=win,nrow=dim(x)[1]+win*2)
  tx<-cbind(ver,rbind(hor,x_1,hor),ver)
  #
  ## Loop over all the pixels
  #
  for (cl in (1+win):(dim(x)[2]+win)) {
    # Update progress bar
    pb$tick()
  for(rw in (1+win):(dim(x)[1]+win)) {
      if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
        out[rw-win,cl-win]<-NA
      } else {
        tw<-summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]))
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if(debugging) {
            message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window^2)
        }
        tw_values<-as.vector(tw)
        out[rw-win,cl-win]<-max(tw_values/sum(tw_values))
      }
    }
  } 
  return(out)
}