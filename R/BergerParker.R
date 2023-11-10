#' Berger-Parker's diversity index
#'
#' Computes Berger-Parker's diversity index on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects. In the latter case, only the first element of the list will be considered.
#' @param window The side of the square moving window, it must be an odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3.
#' @param rasterOut Boolean, if TRUE, output will be in SpatRaster format with \emph{x} as a template.
#' @param np The number of processes (cores) which will be spawned. Default value is 1.
#' @param na.tolerance A numeric value (0.0-1.0) which indicates the proportion of NA values that will be tolerated to calculate Berger-Parker's index in each moving window over \emph{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA, otherwise, Rao's index will be calculated considering the non-NA values. Default values are 1.0 (i.e., no tolerance for NA's).
#' @param cluster.type The type of cluster which will be created. The options are \code{"MPI"} (calls "makeMPIcluster"), \code{"FORK"}, and \code{"SOCK"} (call "makeCluster"). Default type is \code{"SOCK"}.
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For de-bugging only.
#' @details Berger-Parker's index is the relative abundance of the most abundant category (i.e., unique numerical values in the considered numerical matrix). Berger-Parker's index equals the logarithm of the inverse Renyi's index of order infinity, \eqn{log(1/{}^\infty H)} or the inverse of Hill's index of order infinity, \eqn{1/{}^\infty D}.
#' @return A numerical matrix with dimensions as \code{dim(x)}.
#' @references Berger, W.H., Parker, F.L. (1970). Diversity of planktonic foraminifera in deep-sea sediments". Science, 168: 1345-1347.
#' @author Marcantonio Matteo \email{marcantoniomatteo@gmail.com}, Martina Iannacito \email{martina.iannacito@inria.fr}, Duccio Rocchini \email{duccio.rocchini@unibo.it}
#' @note Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
#' apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi
#'
#' Microsoft Windows users may need some additional work to use "MPI", see:\cr
#' \url{https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/}
#'
#' @examples
#' \dontrun{
#' # Minimal example; compute Renyi's index with alpha 1:5 
#' a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
#' berpar <- BergerParker(x=a, window=3)
#' }
#' @export
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach %dopar%

BergerParker <- function(x, window=3, rasterOut=TRUE, np=1, na.tolerance=1.0, cluster.type="SOCK", debugging=FALSE){

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
    rasterm <- x[[1]]
    if( is(rasterm,"SpatRaster") ) {
      rasterm <- matrix(terra::values(rasterm), ncol = ncol(rasterm), nrow = nrow(rasterm), byrow=TRUE)
    }
  }

  #Print user messages
  message("\nObject x check OK: \nBerger-Parker output matrix will be returned.")
    # Derive operational moving window
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of the moving window must be an odd number. Exiting...")
  } 
  if (np == 1){
    outS <- BergerParkerS(rasterm, w, na.tolerance, debugging)
    message(("\nCalculation complete.\n"))
    if(rasterOut==TRUE & class(x)[1]=="SpatRaster") {
      outR <- terra::rast(outS, crs=terra::crs(x),  ext=terra::ext(x))
      return(outR)
    }else{
      return(outS)
    }
  }
  else if (np>1){
  # If more than 1 process
    message("\n##################### Starting parallel calculation #######################")
    if(debugging){cat("#check: Berger-Parker parallel function.")}
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
    outP <- do.call(cbind,BergerParkerP(rasterm, w, na.tolerance, debugging))
    message("\nCalculation complete.\n")
    if(rasterOut==TRUE & class(x)[1]=="SpatRaster") {
      outR <- terra::rast(outP, crs=terra::crs(x),  ext=terra::ext(x))
      return(outR)
    }else{
      return(outP)
    }
  }
}