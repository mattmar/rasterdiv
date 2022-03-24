paRao <- function(x, area=NULL, field=NULL, dist_m="euclidean", window=9, alpha=1, method="classic", rasterOut=TRUE, lambda=0, na.tolerance=1.0, rescale=FALSE, diag=TRUE, simplify=1, np=1, cluster.type="SOCK", debugging=FALSE)
{
    # Define function to check if a number is an integer
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
    # Warning on numeric "simplification"
    message(paste0("Warning: simplify=", simplify,". You're rounding data to ",simplify," decimal place."))
    # Initial checks on type of input data
	if( !(is(x,"matrix") | is(x,"SpatialGridDataFrame") | is(x,"RasterLayer") | is(x,"list")) ) {
		stop("\nNot a valid x object.")
	}
	if( is(x,"SpatialGridDataFrame") ) {
		rasterm <- list(raster(x))
	} else if( is(x,"matrix") | is(x,"RasterLayer") ) {
		rasterm <- list(x)
	} else if( is(x,"list") & method=="multidimension" ) {
		rasterm <- x
	} else if( is(x,"list") & method=="classic" ) {
		stop("If x is a list then method should be *multidimension*?")
	} else stop("Please provide a valid raster input object.")
	if( na.tolerance>1.0 | na.tolerance<0.0 ){
		stop("na.tolerance must be in the [0-1] interval.")
	}
    # Area check
	if( !is.null(area) ) {
		if( !is(area,"SpatialPolygonsDataFrame") ) {stop("area must be SpatialPolygonsDataFrame.")}
		if( !field%in%names(area) ) {stop("field must be a valid variable name of `area`.")}
		if(np>1) {stop("Parallell Area-based Rao's index not implemented yet.")}
		message("Area-based Rao's index")
	}
    # Alpha's check
	if ( any(!is.numeric(alpha)) ){
		stop("alpha must be a numeric vector.")
	}
	if ( any(alpha<0) ){
		stop("alphas must be only positive numbers.")
	}
    # Deal with matrix and RasterLayer in different ways
    # If data are raster layers
	if(is.null(area)){
		if( any(sapply(rasterm, is,"RasterLayer")) ) {
			isfloat <- FALSE
			israst <- TRUE
		# If data are float numbers, transform them to integers. This allows for a shorter computation time on big datasets.
			if( any(sapply(rasterm,function(x) !all(as.matrix(x)==as.integer(as.matrix(x)),na.rm=TRUE))) ){
				message("Input data are float numbers. Converting data to integer matrices...")
				isfloat <- TRUE
				mfactor <- 100^simplify
				rasterm <- lapply(rasterm, function(z) {
					if(rescale) {
						message("Centering and scaling data...")
						z <- as.matrix(raster::scale(z,center=TRUE,scale=TRUE))
					}
					y <- as.matrix(z) * mfactor
					storage.mode(y) <- "integer"
					return(y)
				})
		# If data are integers, just be sure that the storage mode is integer
			}else{
				nr <- sapply(rasterm,nrow); nc <- sapply(rasterm,ncol)
				rasterm <- lapply(rasterm, function(z) 
				{
					if(rescale) {
						message("Centring and scaling data...")
						z <- as.matrix(raster::scale(z,center=TRUE,scale=TRUE))
						mfactor <- 100^simplify
						y <- z * mfactor
						storage.mode(y) <- "integer"
					} else{
						y <- type.convert(matrix(getValues(z),ncol=nc,nrow=nr,byrow=TRUE), as.is= TRUE)
					}
					return(y)
				})
			}
			message("Numerical matrix ready: \nParametric Rao output will be returned")
    # If data are in a matrix or a list
		}else if( any(sapply(rasterm, is,"matrix")) ) {
			isfloat <- FALSE
			israst <- FALSE
		# If data are float numbers, transform them in integer
			if( any(sapply(rasterm, function(x) !all(x == as.integer(x)))) ){
				message("Input data are float numbers. Converting data to integer matrices...")
				isfloat <- TRUE
				mfactor <- 100^simplify
				nr <- sapply(rasterm,nrow); nc <- sapply(rasterm,ncol)
				rasterm <- lapply(rasterm, function(z) {
					if(rescale) {
						message("Centring and scaling data...")
						z <- (z-mean(z))/sd(z)
					}
					y <- z * mfactor
					storage.mode(y) <- "integer"
					return(y)
				})
		# If data are integers, just be sure that the storage mode is integer
			}else{
				rasterm <- lapply(rasterm, function(z) {
					if(rescale) {
						message("Centring and scaling data...")
						z <- (z-mean(z))/sd(z)
						mfactor <- 100^simplify
						y <- z * mfactor
						storage.mode(y) <- "integer"
					}
					type.convert(as.matrix(z), as.is=TRUE)
				})
			}
		}
		message("Numerical matrix ready: \nParametric Rao output will be returned")
	} else ("The class of x is not recognized. Exiting...") 
	# Derive operational moving window
	if( all(window%%2==1) ){
		w <- (window-1)/2
	} else {
		stop("The size of the moving window must be an odd number. Exiting...")
	}
    # Run functions and save outputs
	if( np==1 ) {
		if(method=="classic") {
			if(!is.null(area)) {
				if( debugging ){cat("#check: Inside area clause.")}
				split_layers <- split(area, area[[field]])
				out <- lapply(X=split_layers, function(are){
					lapply(X=alpha, area=are, FUN=paRaoAreaS, rasterm=rasterm[[1]], simplify=simplify)
				})
			} else {
				out <- lapply(X=w, function(win){
					lapply(X=alpha, FUN=paRaoS, rasterm=rasterm[[1]], w=win, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor)
				})}
			} else if(method=="multidimension") {
				out <- lapply(X=w, function(win){
					lapply(X=alpha, FUN=mpaRaoS, x=rasterm, w=win, dist_m=dist_m, na.tolerance=na.tolerance, rescale=rescale, lambda=lambda, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor)
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
					lapply(X=alpha, FUN=paRaoP, rasterm=rasterm[[1]], w=win, dist_m=dist_m, na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor, np = np)
				})
			} else if(method=="multidimension") {
				out <- lapply(X=w, function(win){
					lapply(X=alpha, FUN=mpaRaoP, x=rasterm, w=win, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor, rescale=rescale,np = np)
				})
			}
		}
    # Check if it's an area or moving window based RaoQ
		if( !is.null(area) )
		{
			y <- do.call(rbind.data.frame, lapply(out, function(x) rbind(x)))
			if(nrow(y)>1) y <- as.data.frame(sapply(y,unlist))
			names(y) <- paste("alpha.",alpha, sep="")
			area@data <- cbind.data.frame(area@data,y)
			return(area)
    # Check if the output is a raster or a matrix
		}else{
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
	}