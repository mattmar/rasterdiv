Hill <- function(x, window=3, alpha=1, rasterOut=TRUE, np=1, na.tolerance=1.0, cluster.type="SOCK", debugging=FALSE) {

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
  else if( is(x,"RasterLayer") ) {
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
      if(rasterOut==TRUE & class(x)[1]=="RasterLayer") {
        outR <- lapply(list(outS),raster, template=x)
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
      if(rasterOut==TRUE & class(x)[1]=="RasterLayer") {
        outR <- lapply(outS,raster,template=x)
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
      if(rasterOut==TRUE & class(x)[1]=="RasterLayer") {
        outR <- lapply(list(do.call(cbind,outP)),raster,template=x)
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
      if(rasterOut==TRUE & class(x)[1]=="RasterLayer") {
        outR <- lapply(outP,raster,template=x)
        return(outR)
      }else{
        return(outP)
      }
    }
  }
}