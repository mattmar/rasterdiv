#' Shannon's diversity index (H')
#'
#' Computes Shannon's diversity index (H') on different classes of numeric matrices 
#' using a moving window algorithm.
#'
#' @param x Input data; can be a matrix, a Spatial Grid Data Frame, a SpatRaster, 
#' or a list of these objects. Only the first element of the list is considered if 
#' a list is provided.
#' @param window Numeric; the size of the square moving window. It must be an odd 
#' number greater than 1 to ensure the target pixel is centrally located. Default is 3.
#' @param rasterOut Boolean; if TRUE, the output is a SpatRaster with \code{x} as 
#' the template.
#' @param np Number of processes (cores) to be used; default is 1.
#' @param na.tolerance Numeric between 0.0 and 1.0, indicating the acceptable proportion 
#' of NA values in each moving window over \code{x} for Shannon's index calculation. 
#' If the relative proportion of NAs exceeds na.tolerance, the result is set as NA. 
#' Otherwise, the index is computed using non-NA values. Default is 1.0.
#' @param cluster.type Specifies the type of cluster to create. Options are "MPI" 
#' (utilizing "makeMPIcluster"), "FORK", and "SOCK" (using "makeCluster"). Default 
#' is "SOCK".
#' @param debugging Boolean; if TRUE, additional diagnostic messages are printed, 
#' intended for debugging. Default is FALSE.
#' @return A numerical matrix with the same dimensions as \code{x}, containing 
#' the calculated Shannon's diversity index.
#' @details Shannon's index (H') is calculated as \eqn{H' = -\sum (p_i * log(p_i))}, 
#' where \eqn{p_i} represents the relative abundance of each category within the 
#' numerical matrix, and the sum is computed across all categories.
#' @references Shannon, C.E. (1948). A mathematical theory of communication. Bell 
#' System Technical Journal, 27: 379-423, 623-656.
#' @author Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}, 
#' Martina Iannacito \email{martina.iannacito@@inria.fr}, 
#' Duccio Rocchini \email{duccio.rocchini@@unibo.it}
#' @note For MPI parallel computing, Linux users need to install libopenmpi. On Ubuntu, 
#' you can use the following commands: \code{apt-get update; apt-get upgrade; apt-get install mpi; 
#' apt-get install libopenmpi-dev; apt-get install r-cran-rmpi}. Windows users might 
#' need additional steps for "MPI" usage, more info at: 
#' \url{https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/}
#' @examples
#' \dontrun{
#' # Minimal example: compute Shannon's index
#' a <- matrix(c(10, 10, 10, 20, 20, 20, 20, 30, 30), ncol = 3, nrow = 3)
#' shannon <- Shannon(x = a, window = 3)
#' }
#' @export

Shannon <- function(x, window=3, rasterOut=TRUE, np=1, na.tolerance=1, cluster.type="SOCK", debugging=FALSE){

# Initial checks
  if( !((is(x,"matrix") | is(x,"SpatialGridDataFrame") | is(x,"SpatRaster") | is(x,"list"))) ) {
    stop("\nNot a valid x object. Exiting...")
  }
  else if( is(x,"matrix") ) {
    rasterm <- x
  }
  else if( is(x,"SpatialGridDataFrame") ) {
    rasterm <- terra::rast(x)
  }
  else if( is(x,"SpatRaster")) {
    rasterm <- matrix(terra::values(x), ncol = ncol(x), nrow = nrow(x), byrow=TRUE)
  } 
  else if( is(x,"list") ) {
    message("x is a list, only first element will be taken.")
    if( !((is(x[[1]],"matrix") | is(x[[1]],"SpatialGridDataFrame") | is(x[[1]],"SpatRaster"))) ) {
      stop("The first element of list x is not a valid object. Exiting...")
    }
    rasterm<-x[[1]]
    if( is(rasterm,"SpatRaster") ) {
      rasterm <- matrix(terra::values(rasterm), ncol = ncol(rasterm), nrow = nrow(rasterm), byrow=TRUE)
    }
  }

  #Print user messages
  message("\nObject x check OK: \nShannon output matrix will be returned.")
    # Derive operational moving window
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of the moving window must be an odd number. Exiting...")
  }  
  if (np == 1){
    outS <- ShannonS(rasterm, w, na.tolerance, debugging)
    message("\nCalculation complete.\n")
    if(rasterOut==TRUE & class(x)[[1]]=="SpatRaster") {
      outR <- terra::rast(outS, crs=terra::crs(x),  ext=terra::ext(x))
      return(outR)
    }else{
      return(outS)
    }
  }
  else if (np>1){

    # If more than 1 process
    message("\n##################### Starting parallel calculation #######################")
    if(debugging){cat("#check: Shannon parallel function.")}
    if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
      cls <- makeCluster(np,type=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
    } 
    else if( cluster.type=="MPI" ) {
      cls <- makeCluster(np,outfile="",useXDR=FALSE,methods=FALSE,output="")
    } 
    else {
      message("Wrong definition for cluster.type. Exiting...")
    }
    registerDoParallel(cls)
    # Close clusters on exit
    on.exit(stopCluster(cls))
    # Garbage collection
    gc()
    outP <- do.call(cbind,ShannonP(rasterm, w, na.tolerance, debugging))
    message(("\nCalculation complete.\n"))
    if(rasterOut==TRUE & class(x)[[1]]=="SpatRaster") {
      outR <- terra::rast(outP, crs=terra::crs(x),  ext=terra::ext(x))
      return(outR)
    }else{
      return(outP)
    }
  }
}