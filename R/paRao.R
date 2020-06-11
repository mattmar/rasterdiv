paRao <- function(x, dist_m="euclidean", window=9, alpha=1, method="classic", rasterOut=TRUE, lambda=0, na.tolerance=1.0, rescale=FALSE, diag=TRUE, simplify=2, np=1, cluster.type="SOCK", debugging=FALSE)
{
	#
	## Define function to check if a number is an integer
	#
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
	#
	## Initial checks
	#
    # Initial checks on type of input data
	if( !(is(x,"matrix") | is(x,"SpatialGridDataFrame") | is(x,"RasterLayer") | is(x,"list")) ) {
		stop("\nNot a valid x object.")
	}
	if( is(x,"SpatialGridDataFrame") ) {
		x <- raster(x) # Change x matrix/ces names
	}
	else if( is(x,"matrix") | is(x,"RasterLayer")) {
		rasterm<-x
	} 
	else if( is(x,"list") ) {
		rasterm<-x[[1]]
	}
	if(na.tolerance>1.0 | na.tolerance<0.0){
		stop("na.tolerance must be in the [0-1] interval. Exiting...")
	}
    #Alpha's check
	if ( any(!is.numeric(alpha)) ){
		stop("alpha must be a numeric vector. Exiting...")
	}
	if ( any(alpha<0) ){
		stop("alphas must be only positive numbers. Exiting...")
	}

	#Deal with matrices and RasterLayer in a different way
	#If data are raster layers
	if( method=="classic" & is(x,"RasterLayer") ) {
		isfloat <- FALSE # If data are float numbers, transform them in integer, this may allow for a shorter computation time on big datasets.
		if( !is.wholenumber(rasterm@data@min) | !is.wholenumber(rasterm@data@max) | is.infinite(rasterm@data@min) | !is.wholenumber(median(getValues(rasterm),na.rm=T)) ){
			message("Input data are float numbers. Converting x data in an integer matrix...")
			isfloat<-TRUE
			mfactor<-100^simplify
			rasterm<-getValues(rasterm)*mfactor
			rasterm<-as.integer(rasterm)
			rasterm<-matrix(rasterm,nrow(x),ncol(x),byrow=TRUE)
			gc()
		}else{
			rasterm<-matrix(getValues(rasterm),ncol=ncol(x),nrow=nrow(x),byrow=TRUE)
		}
    #Print user messages
		message("Matrix check OK: \nParametric Rao output matrix will be returned")
    #If data are a matrix or a list
	}else if( method=="classic" & (is(x,"matrix") | is(x,"list")) ) {
		isfloat <- FALSE # If data are float numbers, transform them in integer
		if( !is.integer(rasterm) ){
			message("Input data are float numbers. Converting x in an integer matrix...")
			isfloat<-TRUE
			mfactor<-100^simplify
			rasterm<-as.integer(rasterm*mfactor)
			rasterm<-matrix(rasterm,nrow(x),ncol(x),byrow=TRUE)
			gc()
		}else{
			rasterm<-as.matrix(rasterm)
		}
		message("Matrix check OK: \nParametric Rao output matrix will be returned")
	}else ("The class of x is not recognized. Exiting...") 
	#
	##Derive operational moving window
	#
	if( window%%2==1 ){
		w <- (window-1)/2
	} else {
		stop("The size of the moving window must be an odd number. Exiting...")
	}
	#
	## Run functions and save outputs
	#
	if(np==1) {
		if(method=="classic") {
			out <- pblapply(alpha, paRaoS, rasterm=rasterm, w=w, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat,mfactor=mfactor)
		} else if(method=="multidimensional") {
			out <- pblapply(alpha, mpaRaoS, x=x, rasterm=rasterm, w=w, dist_m=dist_m, na.tolerance=na.tolerance, rescale=rescale, lambda=lambda, diag=diag, debugging=debugging)
		}
    	# Check if the output will be a raster
		if(rasterOut==TRUE & class(x)[[1]]=="RasterLayer") {
			outR <- lapply(out,raster,template=x)
			return(outR)
		}else{
			return(out)
		}
	} else if(np>1) {
		if(method=="multidimensional") {
			stop("Multidimensional paRao not yet implemented as a parallelised function, set 'np=1'. Exiting...")
		} else {
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
				message("Wrong definition for 'cluster.type'. Exiting...")
			}
			doParallel::registerDoParallel(cls)
    		# Close clusters on exit
			on.exit(stopCluster(cls))
    		# Garbage collection
			gc()
			out <- pblapply(alpha, paRaoP, rasterm=rasterm, w=w, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat,mfactor=mfactor)
  			# Check if the output is a raster
			if(rasterOut==TRUE & class(x)[[1]]=="RasterLayer") {
				outR <- lapply(out,raster,template=x)
				return(outR)
			}else{
				return(out)
			}
		}
	}
}