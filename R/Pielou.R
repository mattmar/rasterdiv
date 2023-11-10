#' Pielou's evenness index (E')
#'
#' Computes Pielou's evenness index on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster or a list of these objects. In the latter case, only the first element of the list will be considered.
#' @param window The side of the square moving window, it must be an odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3.
#' @param rasterOut Boolean, if TRUE the output will be in SpatRaster format with \emph{x} as template.
#' @param np The number of processes to use for parallel processing. Defaults to 1.
#' @param na.tolerance A numeric value (0.0-1.0) indicating the proportion of NA values tolerated for calculating Pielou's index in each moving window over \emph{x}. Defaults to 1.0 (i.e., full tolerance).
#' @param cluster.type Specifies the type of cluster to create. Available options are "MPI" (utilizing "makeMPIcluster"), "FORK", and "SOCK" (using "makeCluster"). Default is "SOCK".
#' @param debugging A boolean; if TRUE, additional debugging messages will be printed. Defaults to FALSE.
#'
#' @details Pielou evenness's index (E') is calculated as E' = \deqn{{\sum_{i=1}^{R} p_i \times log(p_i)}\over{log(R)}}, where \emph{R} represents the total number of categories (i.e., unique numerical values) and \emph{p} is the relative abundance of each category. The index represents the ratio between the observed diversity and the maximum possible diversity.
#'
#' @return A numerical matrix with the same dimensions as \code{x}.
#'
#' @references 
#' Pielou, E.C. (1966). The measurement of diversity in different types of biological collections. 
#' Journal of Theoretical Biology, 13: 131-144.
#'
#' @note 
#' Linux users need to install libopenmpi for MPI parallel computing. 
#' For Linux Ubuntu, try the following commands: 
#' apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi
#' Windows users may need to perform additional setup for "MPI". More information is available at:
#' \url{https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/}
#'
#' @seealso \code{\link{Shannon}}
#'
#' @examples
#' # Minimal example; compute Pielou's index  
#' a <- matrix(c(10,10,10,20,20,20,20,30,30), ncol=3, nrow=3)
#' renyi <- Pielou(x=a, window=3)
#'
#' @export

Pielou <- function(x, window=3, rasterOut=TRUE, np=1, na.tolerance=1, cluster.type="SOCK", debugging=FALSE){

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
  message("\nObject x check OK: \nPielou output matrix will be returned.")
    # Derive operational moving window
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of the moving window must be an odd number. Exiting...")
  }  
  if (np == 1){
    outS <- PielouS(rasterm, w, na.tolerance, debugging)
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
    if(debugging){cat("#check: Pielou parallel function.")}
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
    outP <- do.call(cbind,PielouP(rasterm, w, na.tolerance, debugging))
    message(("\nCalculation complete.\n"))
    if(rasterOut==TRUE & class(x)[[1]]=="SpatRaster") {
      outR <- terra::rast(outP, crs=terra::crs(x),  ext=terra::ext(x))
      return(outR)
    }else{
      return(outP)
    }
  }
}