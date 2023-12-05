#' Sequential Pielou's diversity index
#'
#' Computes Pielou's diversity index using a sequential method, particularly useful 
#' for handling large datasets that might not be efficiently processed in a 
#' standard, non-sequential manner.
#'
#' @param x Input raster data, representing the environmental variable(s) 
#' over which the diversity index should be calculated.
#' @param window The size of the half-side of the square moving window used in the 
#' calculation. This determines the scale at which diversity is assessed.
#' @param na.tolerance A numeric value (between 0.0 and 1.0) indicating the 
#' proportion of NA values that are acceptable in each moving window over the 
#' raster data. If the proportion of NA values in a window exceeds this 
#' threshold, the resulting value for that window is set as NA. The default 
#' is 0.0, indicating no tolerance for NA values.
#' @param debugging Boolean flag indicating whether additional console 
#' output should be generated for debugging purposes. Defaults to FALSE.
#'
#' @return A matrix or list of matrices, depending on the input, containing 
#' the calculated Pielou diversity index values. Each cell in the output 
#' matrix represents the diversity index calculated from the corresponding 
#' moving window of the input data.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@gmail.com}, 
#' Martina Iannacito \email{martina.iannacito@inria.fr}, 
#' Duccio Rocchini \email{duccio.rocchini@unibo.it}
#'
#' @seealso \code{\link{Pielou}} for the standard computation of Pielou's 
#' diversity index.
#'
#' @examples
#' \dontrun{
#' # Demonstration of function with hypothetical data
#' # Ensure you replace this with actual raster data
#' demo_raster <- #... (your raster data here)
#' result <- PielouS(x = demo_raster, win = 3, na.tolerance = 0.1, debugging = FALSE)
#' # proceed with analyzing 'result'
#' }

PielouS <- function(x, window = 1, na.tolerance=1, debugging=FALSE){
   # `win` is the operative moving window
   win = window 
   NAwin <- 2*window+1
   message("\n\nProcessing moving Window: ", NAwin)
  # Set a progress bar
  win = 2*window+1
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed",
    total = (dim(x)[2]+NAwin), 
    clear = FALSE, 
    width = 60, 
    force = FALSE)
  # Reshape values
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
      if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA))  < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
        out[rw-win,cl-win]<-NA
        } else {
          tw<-summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]))
          if( "NA's"%in%names(tw) ) {
            tw<-tw[-length(tw)]
          }

          if(debugging) {
            message("\nPielou\nWorking on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",2*win+1)
          }

          tw_values <- as.vector(tw)
          maxS <- log(length(tw))
          p <- tw_values/sum(tw_values)
          p_log <- log(p)
          out[rw-win,cl-win] <- (-(sum(p*p_log)))/maxS

          if( debugging ) {
            message("\ncat: ",paste(tw,collapse=" ")," log S: ",maxS," Pielou: ",out[rw-win,cl-win])
          }

        }
      }
    } 

    return(out)
  }