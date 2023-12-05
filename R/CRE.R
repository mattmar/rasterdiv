#' Cumulative Residual Entropy (CRE) Function
#'
#' @description
#' Computes the Cumulative Residual Entropy (CRE) for spatial raster data. This function 
#' can be used with either a single raster layer or a list of raster layers. It supports 
#' both classic and multidimensional methods for CRE computation. 
#'
#' @param x A matrix, SpatRaster, or a list of SpatRaster objects.
#' @param window The size of the moving window, must be an odd integer.
#' @param method The method for CRE computation, either "classic" or "multidimensional".
#' @param rasterOut Logical, if TRUE, returns a SpatRaster, else returns a matrix.
#' @param rescale Logical, if TRUE, rescales the data before processing.
#' @param na.tolerance A numeric value between 0 and 1, indicating the tolerance level for NA values.
#' @param simplify Integer, the number of decimal places for data rounding in case of float numbers.
#' @param np The number of parallel processes to use.
#' @param cluster.type The type of parallel cluster to use, options are "SOCK", "FORK", or "MPI".
#' @param progBar logical. If TRUE a progress bar is shown.
#' @param debugging Logical, if TRUE, provides additional debugging information during execution.
#'
#' @return Depending on the 'rasterOut' parameter, this function returns either a SpatRaster or a matrix.
#'
#' @examples
#' \dontrun{
#' # For a matrix input:
#' result <- CRE(matrix_data, window=3, method="classic")
#'
#' # For a SpatRaster input:
#' result <- CRE(raster_data, window=3, method="classic", rasterOut=TRUE)
#' }
#'
#' @export

CRE <- function(x, window=3, method="classic", rasterOut=TRUE, rescale=FALSE, na.tolerance=1.0, simplify=2, np=1, cluster.type="SOCK", progBar=TRUE, debugging=FALSE)
{

# Warning for using experimental features
if (method == "multidimension") {
	warning("Multidimension Rao's index is experimental and should be used with caution.")
}

# Warning for data rounding
if (!is.null(simplify)) {
	warning(paste0("Simplify=", simplify, ". Rounding data to ", simplify, " decimal places."))
}

# Validate input data type
if (!(methods::is(x, "matrix") || methods::is(x, "SpatRaster") || methods::is(x, "list") )) {
	stop("\nInvalid input: x must be a matrix, a SpatRaster or a list.")
} 

# Derive operational moving window
if( all(window%%2==1) ){
	w <- (window-1)/2
	} else {
		stop("The size of the moving window must be an odd number. Exiting...")
	}
	NAwin <- 2*window+1

# Processing based on input type and method
if ((methods::is(x, "matrix") || methods::is(x, "SpatRaster")) && method == "classic") {
	rasterm <- list(x)
	} else if ((methods::is(x, "list") || methods::is(x, "SpatRaster")) && method == "multidimension") {
		rasterm <- x
		} else if (methods::is(x, "list") && method != "multidimension") {
			stop("Invalid input: For a list input, method must be set to 'multidimension'.")
			} else {
				stop("Invalid raster input object provided.")
			}

# Validate na.tolerance
if (na.tolerance > 1.0 || na.tolerance < 0.0) {
	stop("na.tolerance must be a value in the [0, 1] interval.")
}

    # Deal with matrices and RasterLayer in a different way
    # If data are raster layers
    if( any(sapply(rasterm, methods::is,"SpatRaster")) ) {
    	isfloat <- FALSE
    	israst <- TRUE
# If data are float numbers, transform them to integers.
if( !any(sapply(rasterm, terra::is.int)) ){
	warning("Input data are float numbers. Converting data to integer matrices.")
	isfloat <- TRUE
	mfactor <- 100^simplify
	rasterm <- lapply(rasterm, function(z) {
		if(rescale) {
			message("Centring and scaling data...")
			z <- terra::as.matrix(terra::scale(z,center=TRUE,scale=TRUE), wide=TRUE)
		}
		y <- terra::as.matrix(z, wide=TRUE) * mfactor
		storage.mode(y) <- "integer"
		return(y)
		})
	} else {
		nr <- sapply(rasterm,nrow); nc <- sapply(rasterm,ncol)
		rasterm <- lapply(rasterm, function(z) 
		{
			if(rescale) {
				message("Centring and scaling data...")
				z <- terra::as.matrix(terra::scale(z,center=TRUE,scale=TRUE), wide=TRUE)
				mfactor <- 100^simplify
				y <- z * mfactor
				storage.mode(y) <- "integer"
				} else{
					y <- utils::type.convert(matrix(terra::values(z),ncol=nc,nrow=nr,byrow=TRUE), as.is= TRUE)
				}
				return(y)
				})
	}
    # If data are a matrix or a list
    } else if( any(sapply(rasterm, methods::is,"matrix")) ) {
    	isfloat <- FALSE
    	israst <- FALSE
# If data are float numbers, transform them in integer
if( !all(sapply(rasterm, function(x) all(apply(x, c(1, 2), is.integer)))) ){
	warning("Input data are float numbers. Converting data to integer matrices...")
	isfloat <- TRUE
	mfactor <- 100^simplify
	rasterm <- lapply(rasterm, function(z) {
		if(rescale & method=="multidimension") {
			message("Centring and scaling data...")
			z <- (z-mean(z))/stats::sd(z)
		}
		y <- round(z * mfactor)
		return(y)
		})
	}else{
		rasterm <- lapply(rasterm, function(z) {
			if(rescale & method=="multidimension") {
				message("Centring and scaling data...")
				z <- (z-mean(z))/stats::sd(z)
				mfactor <- 100^simplify
				y <- round(z * mfactor)
			}
			utils::type.convert(terra::as.matrix(z, wide=TRUE), as.is=TRUE)
			})
	}
	} else ("The class of x is not recognized. Exiting...")

	if( method=="classic" ){
		message("Matrix check OK: \nCumulative Residual Entropy output matrix will be returned")
		}else if( method=="multidimension" ){
			message(("Matrix check OK: \nA matrix with multimension Cumulative Residual Entropy will be returned"))
			}else{
				stop("Matrix check failed: \nNot a valid x | method | distance, please check all these options")
			}

			if(np>1) {
				if(method=="multidimension"){
					message(
						"Multi-core is not supported for multidimensional Rao, proceeding with 1 core...")
					np=1
					}else{
						message("
							##################### Starting parallel calculation #######################")
					}
				}
				pb <- progress::progress_bar$new(
					format = "[:bar] :percent in :elapsed\n",
            	# Total number of ticks is the number of column +NA columns divided the number of processor.
            	total = (dim(x)[2]/np)+(dim(x)[2]*0.05), 
            	clear = FALSE, 
            	width = 60, 
            	force = FALSE)
    #
    ## Preparation of output matrices
    #
    if(np==1) {
    	rasterm<-rasterm[[1]]
    	raoqe<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
    }
    #
    ## If method is classic Rao
    #
    if(method=="classic") {
        #
        ## If classic Cumulative Residual Entropy is parallelized
        #
        if(np>1) {
            #
            ## Reshape values
            #
            rasterm<-rasterm[[1]]
            values<-as.numeric(as.factor(rasterm))
            rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
            #
            ## Add additional columns and rows to match moving window
            #
            hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
            ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
            trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
            rm(hor,ver,rasterm_1,values); gc()
            if(debugging){cat("#check: Cumulative Residual Entropy parallel function.")}
            cls <- openCluster(cluster.type, np, debugging); on.exit(stopCluster(cls)); gc()
            if(debugging){cat("#check: After cluster opening")}
            #
            ## Start the parallelized loop over iter
            #           
            raoqe <- foreach::foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
            	if(debugging) {
            		cat(paste(cl))
            	}
            	# Update progress bar
            	pb$tick()
            	raout <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
            		if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
            			vv<-NA
            			return(vv)
            		} 
            		else {
            			tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            			if( "NA's"%in%names(tw) ) {
            				tw<-tw[-length(tw)]
            			}
            			if( debugging ) {
            				message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window)
            			}
            			tw_labels <- names(tw)
            			tw_values <- as.vector(tw)
                        #if clause to exclude windows with only 1 category
                        if( length(tw_values) <2 ) {
                        	vv<-NA
                        	return(vv)
                        }
                        else {
                        	vv <- .CRE_(tw_values)
                        	return(vv)
                        }
                    }
                    })
            	return(raout)
            } # End classic RaoQ - parallelized
            message(("\n\nCalculation of Cumulative Residual Entropy complete.\n"))
        #
        ## If classic RaoQ is sequential
        #
        } else if(np==1) {
            # Reshape values
            values<-as.numeric(as.factor(rasterm))
            rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
            # Add additional columns and rows for moving window
            hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
            ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
            trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
            # Derive distance matrix
            classes<-levels(as.factor(rasterm))
            # Loop over each pixel
            for (cl in (1+w):(dim(rasterm)[2]+w)) {
            # Update progress bar
            pb$tick()
            for(rw in (1+w):(dim(rasterm)[1]+w)) {
            	if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
            		raoqe[rw-w,cl-w]<-NA
            		} else {
            			tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            			if( "NA's"%in%names(tw) ) {
            				tw<-tw[-length(tw)]
            			}
            			if(debugging) {
            				message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window)
            			}
            			tw_labels <- names(tw)
            			tw_values <- as.vector(tw)
                        #if clause to exclude windows with only 1 category
                        if(length(tw_values) < 2) {
                        	raoqe[rw-w,cl-w]<-NA
                        	} else {
                        		if(isfloat) {
                        			raoqe[rw-w,cl-w]<-.CRE_(tw_values/mfactor)
                        			} else {
                        				raoqe[rw-w,cl-w]<-.CRE_(tw_values)
                        			}
                        		}
                        	} 
                        } 
            } # End of for loop 
            message(("\nCalculation of Cumulative Residual Entropy complete.\n"))
        }
    }  # End classic RaoQ - sequential
    else if( method=="multidimension" ){
    	if(debugging) {
    		message("#check: Into multidimensional clause.")
    	}
        #----------------------------------------------------#
        #
        ## If multimensional Cumulative Residual Entropy
        #
        # Check if there are NAs in the matrices
        if ( methods::is(rasterm,"SpatRaster") ){
        	if(any(sapply(lapply(unlist(x),length),is.na)==TRUE))
        	message("\n Warning: One or more SpatRaster contain NA's which will be treated as 0")
        	} else if ( methods::is(rasterm,"matrix") ){
        		if(any(sapply(x, is.na)==TRUE) ) {
        			message("\n Warning: One or more matrices contain NA's which will be treated as 0")
        		}
        	}
        #
        ## Reshape values
        #
        vls<-lapply(x, function(x) {as.matrix(x)})
        #
        ## Rescale and add additional columns and rows for moving w
        #
        hor<-matrix(NA,ncol=dim(vls[[1]])[2],nrow=w)
        ver<-matrix(NA,ncol=w,nrow=dim(vls[[1]])[1]+w*2)
        if(rescale) {
        	trastersm<-lapply(vls, function(x) {
        		t1 <- terra::scale(terra::rast(cbind(ver,rbind(hor,x,hor),ver)))
        		t2 <- as.matrix(t1)
        		return(t2)
        		})
        	} else {
        		trastersm<-lapply(vls, function(x) {
        			cbind(ver,rbind(hor,x,hor),ver)
        			})
        	}
        	if(debugging) {
        		message("#check: After rescaling in multimensional clause.")
        	}
        #
        ## Loop over all the pixels in the matrices
        #
        if( (ncol(vls[[1]])*nrow(vls[[1]]))> 10000) {
        	message("\n Warning: ",ncol(vls[[1]])*nrow(vls[[1]])*length(vls), " cells to be processed, it may take some time... \n")
        }
        for (cl in (1+w):(dim(vls[[1]])[2]+w)) {
        	    # Update progress bar
        	    pb$tick()
        	    for(rw in (1+w):(dim(vls[[1]])[1]+w)) {
        	    	if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
        	    		raoqe[rw-w,cl-w] <- NA
        	    		} else {
        	    			tw<-lapply(trastersm, function(x) { x[(rw-w):(rw+w),(cl-w):(cl+w)]
        	    				})
                    #
                    ##Vectorize the matrices in the list and calculate
                    #Among matrix pairwase distances
                    lv <- lapply(tw, function(x) {as.vector(t(x))})
                    raoqe[rw-w,cl-w] <- .CRE_(data.frame(lv))
                }
            }
        }
        message("\nCalculation of Multidimensional Cumulative Residual Entropy index complete.\n")
    }
#----------------------------------------------------#

#
## Return multiple outputs; clarity should be improved
#
if(debugging){
	message( "#check: return function." )
}
if( method=="classic" ) {
	if( isfloat & np>1 ) {
		if(rasterOut & class(x)[1]=="SpatRaster") {
			return(terra::rast(do.call(cbind,raoqe)/mfactor, crs=terra::crs(x[[1]]), extent=terra::ext(x)))
			}else{
				return(do.call(cbind,raoqe)/mfactor)
			}
			if(debugging){
				message("#check: return function - classic.")
			}
			} else if( !isfloat & np>1 ) {
				if(rasterOut & class(x)[1]=="SpatRaster") {
					return(terra::rast(raoqe, crs=terra::crs(x[[1]]), extent=terra::ext(x)))
					} else { 
						return(do.call(cbind,raoqe))
					}
					} else {
						if(rasterOut & class(x)[1]=="SpatRaster") {
							return(terra::rast(raoqe, crs=terra::crs(x[[1]]), extent=terra::ext(x)))
							} else {
								return(raoqe)
							} 
						} 
						} else if( method=="multidimension" ) {
							outl <- list(raoqe)
							names(outl) <- c("Multidimension_CRE")
							if(rasterOut & class(x)[1]=="SpatRaster") {
								return(terra::rast(raoqe, crs=terra::crs(x[[1]]), extent=terra::ext(x)))}
								else {
									return(outl)
								}
							}
						}