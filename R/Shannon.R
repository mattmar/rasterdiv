Shannon <- function(x, window=3, np=1, na.tolerance=1, cluster.type="SOCK", debugging=FALSE){

# Initial checks
  if( !((is(x,"matrix") | is(x,"SpatialGridDataFrame") | is(x,"RasterLayer") | is(x,"list"))) ) {
    stop("\nNot a valid x object. Exiting...")
  }
  else if( is(x,"matrix") ) {
    rasterm <- x
  }
  else if( is(x,"SpatialGridDataFrame") ) {
    rasterm <- raster(x)
  }
  else if( is(x,"RasterLayer")) {
    rasterm <- matrix(getValues(x), ncol = ncol(x), nrow = nrow(x), byrow=TRUE)
  } 
  else if( is(x,"list") ) {
    message("x is a list, only first element will be taken.")
    if( !((is(x[[1]],"matrix") | is(x[[1]],"SpatialGridDataFrame") | is(x[[1]],"RasterLayer"))) ) {
      stop("The first element of list x is not a valid object. Exiting...")
    }
    rasterm<-x[[1]]
    if( is(rasterm,"RasterLayer") ) {
      rasterm <- matrix(getValues(rasterm), ncol = ncol(rasterm), nrow = nrow(rasterm), byrow=TRUE)
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
    return(outS)
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
    return(outP)
  }
}