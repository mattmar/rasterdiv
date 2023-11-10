#' Renyi's entropy (H)
#'
#' Computes Renyi's entropy (\deqn{H}) on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data, which may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects. 
#' In the latter case, only the first element of the list will be considered.
#' @param window The side of the square moving window; it must be an odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. 
#' Default value is 3.
#' @param alpha Order of diversity to compute the index. If \code{alpha} is a vector with length greater than 1, then the index will be calculated over \code{x} for each value in the sequence.
#' @param base A numerical value defining the base of the logarithm in Renyi's entropy formula. Default is \code{exp(1)}.
#' @param rasterOut Boolean; if TRUE, the output will be in SpatRaster format with \code{x} as a template.
#' @param np The number of processes (cores) to be spawned. Default is 1.
#' @param na.tolerance A numeric value between 0.0 and 1.0 indicating the proportion of NA values tolerated when calculating Renyi's index in each moving window over \code{x}. 
#' If the relative proportion of NAs is larger than \code{na.tolerance}, the value for that window is set as NA; otherwise, Renyi's index is calculated considering the non-NA values. Default is 1.0 (no tolerance for NAs).
#' @param cluster.type The type of cluster to create. Options are "MPI" (calls "makeMPIcluster"), "FORK", and "SOCK" (calls "makeCluster"). Default is "SOCK".
#' @param debugging Boolean; if TRUE, additional messages will be printed for debugging purposes. Default is FALSE.
#' @return A list of matrices with a length equal to the length of \code{alpha}. If the length of \code{alpha} is 1, then a matrix of dimensions \code{dim(x)} is returned.
#' @details Renyi's entropy (\deqn{H}) is calculated on a numerical matrix as \deqn{H = \frac{1}{(1-q)} \ln(\sum_{i=1}^{R} {p^q}_i)}, where \emph{q} is the considered order of diversity (\code{alpha}), \emph{R} is the total number of categories (i.e., unique numerical values in the considered numerical matrix), and \emph{p} is the relative abundance of each category. 
#' If q=1, \code{Shannon.R} is called to calculate \deqn{H'} instead of the indefinite \deqn{{}^1D}. If \deqn{p > 2*10^9}, then \code{BerkgerParker.R} is called to calculate \deqn{log(1/{}^\infty H)}. Renyi's entropy of low order gives more weight to rarer categories, while higher orders weight more dominant categories.
#' @references Rényi, A. (1970). Probability Theory. North Holland Publishing Company, Amsterdam.
#' @author Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}, Martina Iannacito \email{martina.iannacito@@inria.fr}, Duccio Rocchini \email{duccio.rocchini@@unibo.it}
#' @note For Linux users: to enable MPI parallel computing, install libopenmpi by running the following commands in the terminal: 
#' \code{apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi}.
#' 
#' Windows users might require additional steps to use "MPI". Refer to the detailed guide at: 
#' \url{https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/}
#' @seealso \code{\link{Shannon}}, \code{\link{BergerParker}}
#' @examples
#' # Minimal example: compute Renyi's index with alpha 1:5 
#' a <- matrix(c(10,10,10,20,20,20,20,30,30), ncol=3, nrow=3)
#' renyi <- Renyi(x = a, window = 3, alpha = 1:5)
#' 
#' @export

Renyi <- function(x, window=3, alpha=1, base=exp(1), rasterOut=TRUE, np=1, na.tolerance=1, cluster.type="SOCK", debugging=FALSE){

  # Initial checks
  if( !((is(x,"matrix") | is(x,"SpatialGridDataFrame") | any(is(x,"SpatRaster")) | is(x,"list"))) ) {
    stop("\nNot a valid x object. Exiting...")
  }
  else if( is(x,"matrix") ) {
    rasterm <- x
  }
  else if( is(x,"SpatialGridDataFrame") ) {
    rasterm <- terra::rast(x)
  }
  else if( any(is(x,"SpatRaster"))) {
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
  if ( any(!is.numeric(alpha)) ){
    stop("alpha must be a numeric vector. Exiting...")
  }
  if ( any(alpha<0) ){
    stop("alpha must be only positive numbers. Exiting...")
  }

  # Assign mode according to the length of alpha
  if( length(alpha)==1 ) mode <- "single" else if( length(alpha)>1 ) mode <- "iterative"
  # One single alpha
  if (mode=="single"){
    # If alpha is ~1 then Shannon is calculated
    if ( abs(alpha-1)<.Machine$double.eps ) {
      Shannon <- TRUE
    } else{
      Shannon <- FALSE
    }
    # If alpha approaches positive infinity then Berger-Parker is calculated
    if ( alpha >= .Machine$integer.max ) {
      BergerParker <- TRUE
    } else {
      BergerParker <- FALSE
    }
  }
  # Print messages about output
  message(paste(c("\nObject x check OK: \nRenyi with alpha parameter value in ", alpha," will be returned."),collapse=" "))
  # Derive operational moving window
  if( window%%2==1 ){
    w <- (window-1)/2
  } else {
    stop("The size of the moving window must be an odd number. Exiting...")
  }

  # If one single process
  if (np == 1){
    if(mode == "single") {
      if( BergerParker ) {
        outS <- log(1/BergerParkerS(rasterm, w, na.tolerance, debugging))
      }
      else if( Shannon ) {
        outS <- ShannonS(rasterm, w, na.tolerance, debugging)
      }
      else{
        outS <- RenyiS(rasterm, w, alpha, base, na.tolerance, debugging)
      }
      message(("\nCalculation complete.\n"))
      if(rasterOut==T & any(is(x,"SpatRaster"))) {
        outR <- lapply(list(outS),terra::rast, crs=terra::crs(x), ext=terra::ext(x))
        return(outR)
      }else{
        return(outS)
      }
    }
    else if(mode == "iterative"){
      outS <- list()
      for (ALPHA in alpha){
        message("\nProcessing alpha ",ALPHA,"\n")
        if((abs(ALPHA-1)<.Machine$double.eps)) {
          s <- "Shannon_Renyi_alpha_1"
          outS[[s]] <- ShannonS(rasterm, w, na.tolerance, debugging)
        }
        else if (ALPHA >= .Machine$integer.max) {
          s <- "Berger-Parker"
          outS[[s]] <- log(1/BergerParkerS(rasterm, w, na.tolerance, debugging))
        }
        else{
          s <- paste("Renyi_alpha_",as.character(ALPHA),sep="")
          outS[[s]] <- RenyiS(rasterm, w, ALPHA, base, na.tolerance, debugging)
        }
      }
      message(("\nCalculation complete.\n"))
      if(rasterOut==T & any(is(x,"SpatRaster"))) {
        outR <- lapply(outS, terra::rast, crs=terra::crs(x), ext=terra::ext(x))
        return(outR)
      }else{
        return(outS)
      }
    }
  }
  else if (np>1){

    # If more than 1 process
    message("\n##################### Starting parallel calculation #######################")     
    # Opening the cluster
    if(debugging){cat("#check: Before parallel function.")}
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
    if(mode == "single") {
      if(BergerParker){
        outP <- 1/log(BergerParkerP(rasterm, w, na.tolerance, debugging))
      }
      if( Shannon ) {
        outP <- ShannonP(rasterm, w, na.tolerance, debugging)
      }
      else{
        outP <- RenyiP(rasterm, w, alpha, base, na.tolerance, debugging)
      }
      if(rasterOut==T & any(is(x,"SpatRaster"))) {
        outR <- lapply(list(do.call(cbind,outP)),terra::rast,  terra::crs(x),  terra::res(x))
        return(outR)
      }else{
        return(outP)
      }
    }
    else if(mode == "iterative"){
      outP <- list()
      for (ALPHA in alpha){
        if( (abs(ALPHA-1)<.Machine$double.eps) ) {
          s <- "Shannon_Renyi_alpha_1"
          out <- ShannonP(rasterm, w, na.tolerance, debugging)
          outP[[s]] <- do.call(cbind,out)
        } 
        else if ( ALPHA >= .Machine$integer.max ){
          s <- "Berger-Parker"
          out <- BergerParkerP(rasterm, w, na.tolerance, debugging)
          outP[[s]] <- 1/log(do.call(cbind,out))
        }
        else{
          s <- paste("Renyi_alpha_",as.character(ALPHA),sep="")
          out <- RenyiP(rasterm, w, ALPHA, base, na.tolerance, debugging)
          outP[[s]] <- do.call(cbind,out)
        }
      }
      if(rasterOut==T & any(is(x,"SpatRaster")) ) {
        outR <- lapply(outP,terra::rast, terra::crs(x),  terra::res(x))
        return(outR)
      }else{
        return(outP)
      }
      message("\nCalculation complete.\n")
    }
  }
}