#' Calculate Shannon-Wiener Index on a Matrix
#'
#' @description
#' This function computes Shannon-Wiener Index for each cell of a matrix, 
#' using a parallelized approach and considering a specified moving window.
#'
#' @param x A numeric matrix representing the data on which the index is to abe calculated.
#' @param window The width of the moving window to consider for each cell. 
#'        The actual window size will be `(2 * window + 1) x (2 * window + 1)`. Default is 1.
#' @param na.tolerance The tolerance level for missing data within the moving window. 
#'        A window will be processed only if the proportion of non-missing data is above this threshold. 
#'        Value should be between 0 and 1. Default is 1.
#' @param debugging Boolean flag to enable or disable debugging messages. Default is FALSE.
#' @param np Number of processes for parallel computation.
#'
#' @return A matrix of the same dimensions as `x`, where each cell contains the 
#'         Shannon-Wiener Index calculated for the window around the cell.
#'
#' @examples
#' data <- matrix(runif(100), nrow = 10)
#' shannon_index <- ShannonP(data, window = 1, np = 1 )
#'
#' @export

ShannonP <- function(x, window = 1, na.tolerance=1, debugging=FALSE, np =1 ){
  # `win` is the operative moving window
  win = window 
  NAwin <- 2*window+1
  message("\n\nProcessing moving Window: ", NAwin)
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
  values <- as.numeric( as.factor(x) )
  x_1 <- matrix(data = values, nrow = nrow(x), ncol = ncol(x))
  #
  ## Add additional columns and rows to match moving window
  #
  hor <- matrix(NA, ncol = ncol(x), nrow = win)
  ver <- matrix(NA, ncol = win, nrow = nrow(x)+ win * 2)
  tx <- cbind(ver, rbind(hor,x_1,hor), ver)
  rm(hor, ver, x_1, values); gc()
  ShannonOP <- foreach::foreach(cl=(1+win):(ncol(x)+win),.verbose = F) %dopar% {
    # Update progress bar
    pb$tick()

    ShannonOut <- sapply((1+win):(nrow(x)+win), function(rw) {
      if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
        vv <- NA
        return(vv)
      } 
      else {
        tw <- summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]),maxsum=10000)
        if( "NA's"%in%names(tw) ) {
          tw <- tw[-length(tw)]
        }
        if( debugging ) {
          message("Shannon - parallelized\nWorking on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",NAwin)
        }
        tw_labels <- names(tw)
        tw_values <- as.vector(tw)
        p <- tw_values/sum(tw_values)
        vv <- (-(sum(p*log(p))))
        return(vv)
      }
    })
    return(ShannonOut)
  }
  message("\n\n Parallel calculation of Shannon's index complete.\n")
  return(matrix(unlist(ShannonOP), ncol = ncol(x), nrow = nrow(x), byrow=FALSE))
}