paRao <- function(x, dist_m="euclidean", window=9, alpha=1, method="classic", rasterOut=TRUE, lambda=0, na.tolerance=1.0, rescale=FALSE, diag=TRUE, simplify=1, np=1, cluster.type="SOCK", debugging=FALSE)
{
    # Define function to check if a number is an integer
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
    # Initial checks on type of input data
	if( !(is(x,"matrix") | is(x,"SpatialGridDataFrame") | is(x,"RasterLayer") | is(x,"list")) ) {
		stop("\nNot a valid x object.")
	}
	if( is(x,"SpatialGridDataFrame") ) {
		rasterm <- list(raster(x))
	}
	else if( is(x,"matrix") | is(x,"RasterLayer") ) {
		rasterm <- list(x)
	} 
	else if( is(x,"list") & method=="multidimension" ) {
		rasterm <- x
	}
	else if( is(x,"list") & method=="classic" ) {
		stop("If x is a list then method must be *multidimensional*")
	}
	if( na.tolerance>1.0 | na.tolerance<0.0 ){
		stop("na.tolerance must be in the [0-1] interval. Exiting...")
	}
    # Alpha's check
	if ( any(!is.numeric(alpha)) ){
		stop("alpha must be a numeric vector. Exiting...")
	}
	if ( any(alpha<0) ){
		stop("alphas must be only positive numbers. Exiting...")
	}
    # Deal with matrix and RasterLayer in different ways
    # If data are raster layers
	if( any(sapply(rasterm, is,"RasterLayer")) ) {
		isfloat <- FALSE
		israst <- TRUE
		# If data are float numbers, transform them to integers. This allows for a shorter computation time on big datasets.
		if( any(sapply(rasterm,function(x) grepl("FLT", dataType(x)))) ){
			message("Input data are float numbers. Converting data to integer matrices...")
			isfloat <- TRUE
			nr <- sapply(rasterm,nrow); nc <- sapply(rasterm,ncol)
			mfactor <- 100^simplify
			rasterm <- lapply(rasterm, function(x) {
				y <- type.convert(raster::as.matrix((x * mfactor)))
				return(y)
			})
			gc()
		}else{
			nr <- sapply(x,nrow); nc <- sapply(x,ncol)
			rasterm <- lapply(rasterm, function(x) type.convert(matrix(getValues(x),ncol=nr,nrow=nc,byrow=TRUE)))
		}
		message("Numerical matrix ready: \nParametric Rao output will be returned")
    # If data are a matrix or a list
	}else if( any(sapply(rasterm, is,"matrix")) ) {
		isfloat <- FALSE # If data are float numbers, transform them in integer
		israst <- FALSE
		if( any(sapply(rasterm, function(x) !is.integer(x))) ){
			message("Input data are float numbers. Converting data to integer matrices...")
			isfloat <- TRUE
			mfactor <- 100^simplify
			nr <- sapply(rasterm,nrow); nc <- sapply(rasterm,ncol)
			rasterm <- lapply(rasterm, function(x) {
				y <- type.convert(x * mfactor)
			})
		}else{
			rasterm <- lapply(rasterm, function(x) type.convert(as.matrix(x)))
		}
		message("Numerical matrix ready: \nParametric Rao output will be returned")
	}else ("The class of x is not recognized. Exiting...") 
	# Derive operational moving window
	if( all(window%%2==1) ){
		w <- (window-1)/2
	} else {
		stop("The size of the moving window must be an odd number. Exiting...")
	}
    # Run functions and save outputs
	if( np==1 ) {
		if(method=="classic") {
			out <- lapply(X=w, function(win){
				lapply(X=alpha, FUN=paRaoS, rasterm=rasterm[[1]], w=win, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat,mfactor=mfactor)
			})
		} else if(method=="multidimension") {
			out <- lapply(X=w, function(win){
				lapply(X=alpha, FUN=mpaRaoS, x=rasterm, w=win, dist_m=dist_m, na.tolerance=na.tolerance, rescale=rescale, lambda=lambda, diag=diag, debugging=debugging,isfloat=isfloat, mfactor=mfactor)
			})
		}
	} else if( np>1 ) {
		message("\n##################### Starting parallel calculation #######################")
		if( debugging ){cat("#check: Before parallel function.")}     
        # Opening the cluster
		if( cluster.type=="SOCK" || cluster.type=="FORK" ) {
			cls <- invisible(parallel::makeCluster(np,type=cluster.type, outfile= ' '))
		} 
		else if( cluster.type=="MPI" ) {
			cls <- invisible(makeCluster(np,outfile="",useXDR=FALSE,methods=FALSE,output=""))
		} 
		else {
			message("Wrong definition for 'cluster.type'. Exiting...")
		}
		doParallel::registerDoParallel(cls)
        # Close clusters on exit
		on.exit(stopCluster(cls))
		gc()
		if( method=="classic" ) {
			out <- lapply(X=w, function(win){
				lapply(X=alpha, FUN=paRaoP, rasterm=rasterm[[1]], w=win, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor, np = np)
			})
		} else if(method=="multidimension") {
			out <- lapply(X=w, function(win){
				lapply(X=alpha, FUN=mpaRaoP, x=rasterm, w=win, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor, rescale=rescale,np = np)
			})
		}
	}
    # Check if the output is a raster or a matrix
	if( rasterOut & israst ) {
		outR <- lapply(out, function(insm) {
			y <- lapply(insm,raster,template=x[[1]])
			names(y) <- paste("alpha.",alpha, sep="")
			return(y)
		})
		names(outR) <- paste("window.",window, sep="")
		return(outR)
	}else{
		outM <- lapply(out, function(insm) {
			names(insm) <- paste("alpha.",alpha, sep="")
			return(insm)
		})
		names(outM) <- paste("window.",window, sep="")
		return(outM)
	}
}
