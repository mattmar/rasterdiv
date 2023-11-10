#' Hill's index of diversity - Hill numbers (D)
#'
#' Computes Hill's index of diversity (Hill numbers) on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects. In the latter case, only the first element of the list will be considered.
#' @param window The side of the square moving window. It must be an odd numeric value greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3.
#' @param alpha Order of the Hill number to compute the index. If \code{alpha} is a vector with length greater than 1, then the index will be calculated over \code{x} for each value in the sequence.
#' @param rasterOut Boolean; if TRUE, the output will be in SpatRaster format with \code{x} as the template.
#' @param np The number of processes (cores) which will be spawned. Default value is 1.
#' @param na.tolerance A numeric value between 0.0 and 1.0, which indicates the proportion of NA values that will be tolerated to calculate Hill's index in each moving window over \code{x}. If the relative proportion of NA's in a moving window is bigger than na.tolerance, then the value of the window will be set as NA; otherwise, Hill's index will be calculated considering the non-NA values. Default value is 1.0 (i.e., full tolerance for NA's).
#' @param cluster.type The type of cluster which will be created. Options are "MPI" (calls "makeMPIcluster"), "FORK," and "SOCK" (call "makeCluster"). Default type is "SOCK".
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed for debugging purposes.
#'
#' @details
#' Hill numbers (\eqn{{}^qD}) are calculated on numerical matrices as \eqn{{}^qD = (\sum_{i=1}^{R} {p^q}_i)^{1/(1-q)}}, where \emph{q} is the order of the Hill number, \emph{R} is the total number of categories (i.e., unique numerical values in a numerical matrix), and \emph{p} is the relative abundance of each category. When q=1, Shannon.R is called to calculate \eqn{exp(H^1)} instead of the indefinite \eqn{{}^1D}. If \eqn{q > 2*10^9}, BergerParker.R is called to calculate \eqn{1/{{}^\infty D}}. Hill numbers of low order weight more rare categories, whereas Hill numbers of higher order weight more dominant categories.
#'
#' @return A list of matrices of dimension \code{dim(x)} with length equal to the length of \code{alpha}.
#'
#' @references
#' Hill, M.O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology 54, 427-432.
#'
#' @note
#' Linux users need to install libopenmpi for MPI parallel computing. Linux Ubuntu users may try:
#' \code{apt-get update; apt-get upgrade; apt-get install mpi; apt-get install libopenmpi-dev; apt-get install r-cran-rmpi}
#' 
#' Microsoft Windows users may need some additional work to use "MPI". For more details, see:
#' \url{https://bioinfomagician.wordpress.com/2013/11/18/installing-rmpi-mpi-for-r-on-mac-and-windows/}
#'
#' @seealso
#' \code{\link{BergerParker}}, \code{\link{Shannon}}
#'
#' @examples
#' # Minimal example; compute Hill's index with alpha 1:5 
#' a <- matrix(c(10,10,10,20,20,20,20,30,30),ncol=3,nrow=3)
#' hill <- Hill(x=a,window=3,alpha=1:5)
#'
#' @export

Hill <- function(x, window=3, alpha=1, rasterOut=TRUE, np=1, na.tolerance=1.0, cluster.type="SOCK", debugging=FALSE) {

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
  else if( is(x,"SpatRaster") ) {
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
  if(length(alpha)==1) mode <- "single" else if(length(alpha)>1) mode <- "iterative"
  # One single alpha
  if (mode=="single"){
    # If alpha is ~1 then Shannon is calculated
    if (abs(alpha-1)<.Machine$double.eps) {
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
  message(paste(c("\nObject x check OK: \nHill with alpha parameter value in ", alpha," will be returned."),collapse=" "))
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
        outS <- 1/BergerParkerS(rasterm, w, na.tolerance, debugging)
      }
      else if( Shannon ) {
        outS <- exp(ShannonS(rasterm, w, na.tolerance, debugging))
      }
      else{
        outS <- HillS(rasterm, w, alpha, na.tolerance, debugging)
      }
      message(("\nCalculation complete.\n"))
      # Check if the output will be a raster
      if(rasterOut==TRUE & any(is(x,"SpatRaster"))) {
        outR <- lapply(list(outS),terra::rast,  crs=terra::crs(x),  ext=terra::ext(x))
        return(outR)
        }else{
          return(outS)
        }
      }
      else if(mode == "iterative") {
        outS <- list()
        for (ALPHA in alpha){
          message("\n\nProcessing alpha ",ALPHA)
          if((abs(ALPHA-1)<.Machine$double.eps)) {
            s <- "Shannon_Hill_alpha_1"
            outS[[s]] <- exp(ShannonS(rasterm, w, na.tolerance, debugging))
          }
          else if (ALPHA >= .Machine$integer.max) {
            s <- "Berger-Parker"
            outS[[s]] <- 1/BergerParkerS(rasterm, w, na.tolerance, debugging)
          }
          else{
            s<-paste("Hill_alpha_",as.character(ALPHA),sep="")
            outS[[s]] <- HillS(rasterm, w, ALPHA, na.tolerance, debugging)
          }
        }
        message(("\n\nCalculation complete.\n"))
        if(rasterOut==TRUE & any(is(x,"SpatRaster"))) {
          outR <- lapply(outS,terra::rast, crs=terra::crs(x),  ext=terra::ext(x))
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
    doParallel::registerDoParallel(cls)
    # Close clusters on exit
    on.exit(stopCluster(cls))
    # Garbage collection
    gc()
    if(mode == "single") {
      if( BergerParker ){
        outP <- BergerParkerP(rasterm, w, na.tolerance, debugging)
      }
      else if( Shannon ){
        outP <- lapply(ShannonP(rasterm, w, na.tolerance, debugging),exp)
      }
      else{
        outP <- HillP(rasterm, w, alpha, na.tolerance, debugging)
      }
      message("\nCalculation complete.\n")
      if(rasterOut==TRUE & any(is(x,"SpatRaster"))) {
        outR <- lapply(list(do.call(cbind,outP)),terra::rast, crs=terra::crs(x),  ext=terra::ext(x))
        return(outR)
        }else{
          return(outP)
        }
      }
      else if(mode == "iterative"){
        outP <- list()
        for (ALPHA in alpha){
          if( (abs(ALPHA-1)<.Machine$double.eps) ) {
            s <- "Shannon_Hill_alpha_1"
            out <- ShannonP(rasterm, w, na.tolerance, debugging)
            outP[[s]] <- exp(do.call(cbind,out))
          } 
          else if ( ALPHA >= .Machine$integer.max ){
            s <- "Berger-Parker"
            out <- BergerParkerP(rasterm, w, na.tolerance, debugging)
            outP[[s]] <- do.call(cbind,out)
          }
          else{
            s <- paste("Hill_alpha_",as.character(ALPHA),sep="")
            out <- HillP(rasterm, w, ALPHA, na.tolerance, debugging)
            outP[[s]] <- do.call(cbind,out)
          }
        }
        message("\nCalculation complete.\n")
        if(rasterOut==TRUE & any(is(x,"SpatRaster"))) {
          if(!is.list(outP)) {outP <- list(outP)}
          outR <- lapply(outP,terra::rast, crs=terra::crs(x),  ext=terra::ext(x))
          return(outR)
          }else{
            return(outP)
          }
        }
      }
    }