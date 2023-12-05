#' Calculate Berger-Parker Index on a Matrix
#'
#' @description
#' This function computes Berger-Parker Index for each cell of a matrix, 
#' using a parallelized approach and considering a specified moving window.
#'
#' @param x A numeric matrix representing the data on which the index is to be calculated.
#' @param window The width of the moving window to consider for each cell. 
#'        The actual window size will be `(2 * window + 1) x (2 * window + 1)`. Default is 1.
#' @param na.tolerance The tolerance level for missing data within the moving window. 
#'        A window will be processed only if the proportion of non-missing data is above this threshold. 
#'        Value should be between 0 and 1. Default is 1.
#' @param debugging Boolean flag to enable or disable debugging messages. Default is FALSE.
#' @param np The number of processes (cores) which will be spawned. Default value is 2.
#'
#' @return A matrix of the same dimensions as `x`, where each cell contains the 
#'         Berger-Parker Index calculated for the window around the cell.
#'
#' @examples
#' \dontrun{
#' data <- matrix(runif(100), nrow = 10)
#' bp_index <- BergerParkerP(data, window = 1, np=2)
#' }
#'
#' @export
BergerParkerP<-function(x, window = 1, na.tolerance=1, debugging=FALSE, np=1){
  # `win` is the operative moving window
  win = window 
  NAwin <- 2*window+1
  message("\n\nProcessing moving Window: ", NAwin)
  # Set a progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed\n",
    # Total number of ticks is the number of column +NA columns divided the number of processor.
    total = (dim(x)[2]/np)+5, 
    clear = FALSE, 
    width = 60, 
    force = FALSE)
  #
  ## Reshape values
  #
  values<-as.numeric(as.factor(x))
  x_1<-matrix(data=values,nrow=dim(x)[1],ncol=dim(x)[2])
  #
  ## Add additional xcolumns and rows to match moving window
  #
  hor<-matrix(NA,ncol=dim(x)[2],nrow=win)
  ver<-matrix(NA,ncol=win,nrow=dim(x)[1]+win*2)
  tx<-cbind(ver,rbind(hor,x_1,hor),ver)
  rm(hor,ver,x_1,values); gc()
  #
  ## Start the parallelized loop over iter
  #
  BergerParkerOP <- foreach::foreach(cl=(1+win):(dim(x)[2]+win),.verbose = F) %dopar% {
    if(debugging) {
      cat(paste(cl))
    }
    BergerParkerOut <- sapply((1+win):(dim(x)[1]+win), function(rw) {
      # Update progress bar
      pb$tick()
      #check number of NA's against tolerad number
      if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
        vv<-NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw<-tw[-length(tw)]
        }
        if( debugging ) {
          message("Berger-Parker - parallelised\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=", NAwin)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        vv <- max(tw_values/sum(tw_values))
      }
    })
    return(BergerParkerOut)
  }
  message("\n\n Parallel calculation of Berger Parker's index complete.\n")
  return(matrix(unlist(BergerParkerOP), ncol = ncol(x), nrow = nrow(x), byrow=FALSE))
}
