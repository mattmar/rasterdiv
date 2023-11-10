#' Parametric Rao's index of quadratic entropy (Q)
#'
#' It computes the parametric version of Rao's index of quadratic entropy (Q) on different classes of numeric matrices using a moving window algorithm.
#'
#' @param x Input data may be a matrix, a Spatial Grid Data Frame, a SpatRaster, or a list of these objects.
#' @param area Input vector area layer for area-based calculation.
#' @param field Column name of the vector area layer to use to calculate the index.
#' @param dist_m Define the type of distance to be calculated between numerical categories. `dist_m` can be a character string which defines the name of the distance to derive such as "euclidean". The distance names allowed are the same as for \code{proxy::dist}. Alternatively, `dist_m` can be a function which calculates a user-defined distance, (i.e., \code{function(x,y) {return(cos(y-x)-sin(y-x))}}) or a matrix of distances. If `method="multidimension"` then only "euclidean", "manhattan", "canberra", "minkowski" and "mahalanobis" can be used. Default value is "euclidean". If `dist_m` is a matrix, then the function will assume that the matrix contains the distances.
#' @param window The side of the square moving window, it must be a vector of odd numeric values greater than 1 to ensure that the target pixel is in the centre of the moving window. Default value is 3. `window` can be a vector with length greater than 1, in this case, Rao's index will be calculated over `x` for each value in the vector.
#' @param alpha Weight for the distance matrix. If `alpha = 0`, distances will be averaged with a geometric average, if `alpha=1` with an arithmetic mean, if `alpha = 2` with a quadratic mean, `alpha = 3` with a cubic mean, and so on. if `alpha` tends to infinite (i.e., higher than the maximum integer allowed in R) or `alpha=Inf`, then the maximum distance will be taken. `alpha` can be a vector with length greater than 1, in this case, Rao's index will be calculated over `x` for each value in the vector.
#' @param method Currently, there are two ways to calculate the parametric version of Rao's index. If `method="classic"`, then the normal parametric Rao's index will be calculated on a single matrix. If `method="multidimension"` (experimental!), a list of matrices must be provided as input. In the latter case, the overall distance matrix will be calculated in a multi- or hyper-dimensional system by using the distance measure defined through the function argument `dist_m`. Each pairwise distance is then multiplied by the inverse of the squared number of pixels in the considered moving window, and the Rao's Q is finally derived by applying a summation. Default value is `"classic"`.
#' @param rasterOut Boolean, if TRUE the output will be a SpatRaster object with `x` as a template.
#' @param lambda The value of the lambda of Minkowski's distance. Considered only if `dist_m = "minkowski"` and `method="multidimension"`. Default value is 0.
#' @param na.tolerance Numeric value (0.0-1.0) which indicates the proportion of NA values that will be tolerated to calculate Rao's index in each moving window over `x`. If the relative proportion of NA's in a moving window is bigger than `na.tolerance`, then the value of the window will be set as NA, otherwise Rao's index will be calculated considering the non-NA values. Default value is 1.0.
#' @param rescale Boolean. Considered only if `method="multidimension"`. If TRUE, each element of `x` is rescaled and centred.
#' @param diag Boolean. If TRUE then the diagonal of the distance matrix is filled with 0's, otherwise with NA's. If `diag=TRUE` and `alpha=0`, the output matrix will inexorably be 0's.
#' @param simplify Number of decimal places to be retained to calculate distances in Rao's index. Default `simplify=0`.
#' @param np The number of processes (cores) which will be spawned. Default value is 2.
#' @param cluster.type The type of cluster which will be created. The options are `"MPI"` (which calls "makeMPIcluster"), `"FORK"`, and `"SOCK"` (which call "makeCluster"). Default type is `"SOCK"`.
#' @param debugging A boolean variable set to FALSE by default. If TRUE, additional messages will be printed. For debugging only.
#' @return A list of matrices of dimension `dim(x)` with length equal to the length of `alpha`. If `rasterOut=TRUE` and `x` is a SpatRaster, then the output is a list of SpatRaster objects.
#' @details The parametric Rao's Index (Q) is an extension of Rao's Index which considers a generalized mean between distances. The general formula for the parametric Rao's index is Q_a = \deqn{Q = \sum_{i, j} p_i p_j d_{ij}^{\alpha}}. Where `N` is the number of numerical categories, `i` and `j` are pair of numerical categories in the same moving window, and `alpha` is a weight given to distances. In the "multidimension" Rao's index, first the distances among categories are calculated considering more than one feature, and then the overall Rao's Q is derived by using these distances.
#' @references 
#' Rao, C. R. (1982). Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology, 21(1), 24-43. 
#' 
#' @examples
#' \dontrun{
#' # loading data
#' data(volcano)
#' r <- terra::rast(volcano)
#' 
#' # we want to compute Rao's index on this data using a 3x3 window
#' res <- paRao(x = r, window = 3, alpha = 2, method = "classic")
#' terra::plot(res[[1]][[1]])
#' }
#'
#' @export

paRao <- function(x, area=NULL, field=NULL, dist_m="euclidean", window=9, alpha=1, method="classic", rasterOut=TRUE, lambda=0, na.tolerance=1.0, rescale=FALSE, diag=TRUE, simplify=0, np=1, cluster.type="SOCK", debugging=FALSE) {

# Define function to check if a number is an integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
# Warning on numeric "simplification" and multimensional Rao's
message(paste0("Warning: simplify=", simplify,". You're rounding data to ",simplify," decimal place."))
if( method=="multidimension" ) {warning("Multidimension Rao's index is experimental. Use with caution!")}
# Initial checks on type of input data
if( !(is(x,"matrix") | is(x,"SpatialGridDataFrame") | is (x,"SpatRaster") | is(x,"list") | is(x,"RasterStack")) ) { 
	stop("\nNot a valid xobject.") 
} 
if( is(x,"SpatialGridDataFrame") ) {
	rasterm <- list(terra::rast(x)) 
	} else if( is(x,"matrix") | is(x,"SpatRaster") & method=="classic" ){ 
		rasterm <- list(x) 
		} else if( is(x,"list") | is(x,"SpatRaster") | is(x,"RasterStack") & method=="multidimension" ) {
			rasterm <- x 
			} else if( is(x,"list") & method=="classic" ) { 
				stop ("If x is a list then method should be *multidimension*?") 
				} else stop ("Please provide a valid raster input object.") 
				if( na.tolerance>1.0 | na.tolerance<0.0 ) { 
					stop("na.tolerance must be in the [0-1] interval.") 
				}
# Area check
if( !is.null(area) ) {
	if( !is(area,"SpatialPolygonsDataFrame") & !is(area, "SpatVector") ) {stop("area must be SpatialPolygonsDataFrame|SpatVector.")}
	if( !field%in%names(area) ) {stop("field must be a valid variable name of `area`.")}
	if( np>1 ) {stop("Parallell Area-based Rao's index not implemented yet.")}
	message("Area-based Rao's index")
}
# Alpha's check
if ( any(!is.numeric(alpha)) ){
	stop("alpha must be a numeric vector.")
}
if ( any( alpha<0 ) ){
	stop("alphas must be only positive numbers.")
}
# Deal with matrix and SpatRaster in different ways
# If data are raster layers
if( is.null(area) ){
	if( any(sapply(rasterm, is,"SpatRaster")) ) {
		isfloat <- FALSE
		israst <- TRUE
# If data are float numbers, transform them to integers.
if( !any(sapply(rasterm, terra::is.int)) ){
	message("Input data are float numbers. Converting data to integer matrices...")
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
# If data are integers, just be sure that the storage mode is integer
}else{
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
# If data are in a matrix or a list
}else if( any(sapply(rasterm, is,"matrix")) ) {
	isfloat <- FALSE
	israst <- FALSE
# If data are float numbers, transform them in integer
if( !all(sapply(rasterm, function(x) all(apply(x, c(1, 2), is.integer)))) ){
	message("Input data are float numbers. Converting data to integer matrices...")
	isfloat <- TRUE
	mfactor <- 100^simplify
	rasterm <- lapply(rasterm, function(z) {
		if(rescale) {
			message("Centring and scaling data...")
			z <- (z-mean(z))/stats::sd(z)
		}
		y <- round(z * mfactor)
		return(y)
		})
# If data are integers, just be sure that the storage mode is integer
}else{
	rasterm <- lapply(rasterm, function(z) {
		if(rescale) {
			message("Centring and scaling data...")
			z <- (z-mean(z))/stats::sd(z)
			mfactor <- 100^simplify
			y <- round(z * mfactor)
		}
		utils::type.convert(terra::as.matrix(z, wide=TRUE), as.is=TRUE)
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
	if( method=="classic" ) {
		if( !is.null(area) ) {
			if( debugging ){ cat("#check: Inside classic area clause.") }
			split_layers <- terra::split(area, field)
			out <- lapply(X=split_layers, function(are){
				lapply(X=alpha, area=are, FUN=paRaoAreaS, rasterm=rasterm[[1]], simplify=simplify)
				})
			} else {
				out <- lapply(X=w, function(win){
					lapply(X=alpha, FUN=paRaoS, rasterm=rasterm[[1]], w=win, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor)
					})
			}
			} else if( method=="multidimension" ) {
				if( !is.null(area) ) {
					if( debugging ){ cat("#check: Inside multi area clause.") }
					split_layers <- terra::split(area, field)
					out <- lapply(X=split_layers, function(are){
						lapply(X=alpha, area=are, FUN=mpaRaoAreaS, dist_m=dist_m, rasterm=rasterm, simplify=simplify)
						})
					} else {
						out <- lapply(X=w, function(win){
							lapply(X=alpha, FUN=mpaRaoS, x=rasterm, w=win, dist_m=dist_m, na.tolerance=na.tolerance, rescale=rescale, lambda=lambda, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor)
							})
					}
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
		lapply(X=alpha, FUN=paRaoP, rasterm=rasterm[[1]], w=win, dist_m=dist_m, na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor)
		})
	} else if(method=="multidimension") {
		out <- lapply(X=w, function(win){
			lapply(X=alpha, FUN=mpaRaoP, x=rasterm, w=win, dist_m=dist_m,na.tolerance=na.tolerance, diag=diag, debugging=debugging, isfloat=isfloat, mfactor=mfactor, rescale=rescale)
			})
	}
}
# Check if it's an area or moving window based RaoQ
if( !is.null(area) ) {
	y <- do.call(rbind.data.frame, lapply(out, function(x) rbind(x)))
	if(nrow(y)>1) y <- as.data.frame(sapply(y,unlist))
	names(y) <- paste("alpha.",alpha, sep="")
	terra::values(area) <- cbind.data.frame(area,y)
	return(area)
# Check if the output is either a raster or a matrix
}else{
	if( rasterOut & israst ) {
		outR <- lapply(out, function(insm) {
			if(method=="multidimension"){
				y <- lapply(insm, terra::rast, crs=terra::crs(x[[1]]), ext=terra::ext(x[[1]]))
				} else{
					y <- lapply(insm, terra::rast, crs=terra::crs(x), ext=terra::ext(x))
				}
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