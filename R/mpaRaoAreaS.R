#' Area-Based Sequential Parametric Rao's index of quadratic entropy (Q)
#'
#' Calculates an area-based sequential version of the parametric Rao's index of quadratic entropy (Q). 
#' This function is designed for situations where the diversity index needs to consider geographical 
#' areas and works with raster data representing the distribution of species or other measures.
#'
#' @param rasterm Raster; the input raster data representing variables across a geographic space.
#' @param area Numeric; the input vector data representing the areas of interest.
#' @param alpha Numeric; alpha value for order of diversity in Hill's Index.
#' @param simplify Numeric; the parameter that determines the rounding off of the calculations.
#' @param dist_m Character; type of distance metric used (e.g., "euclidean", "manhattan", etc.).
#' @param rescale Logical; whether to scale and centre the values in each element of the raster data.
#' @param lambda Numeric; lambda parameter for Minkowski distance calculation.
#' @param window Numeric; defines the size of the moving window for the analysis.
#' @return A vector similar to the input, with additional columns representing Rao's index values for 
#' each area.
#' @seealso \code{\link{paRao}} for a related function dealing with the parallel computation of Rao's index.
#' @author Matteo Marcantonio \email{marcantoniomatteo@@gmail.com},
#' Duccio Rocchini \email{duccio.rocchini@@unibo.it},
#' Michele Torresani \email{michele.torresani@@unibo.it}

mpaRaoAreaS <- function(rasterm, area, alpha, simplify, dist_m, rescale, lambda, window) {
# Some initial housekeeping
if(alpha<=0) {
    stop("Alpha<=0 not yet implemented for areas.")
}
mfactor <- 100^simplify
crop1 <- terra::crop(rasterm, area)
crop1dt <- terra::as.matrix(crop1)*mfactor
storage.mode(crop1dt) <- "integer"

# Check for only 1 value in the matrix
if( any(is.na(crop1dt)) & all(apply(crop1dt, 2, function(a) length(unique(a))<=2)) ) {
    mpaRaoOareaS <- NA
    } else {
        # Evaluate Rao's method given alpha
        classes <- levels(as.factor(crop1dt))
        window <- nrow(crop1dt) #Temporary patch?
        if( alpha >= .Machine$integer.max | is.infinite(alpha) ) {
            alphameth <- "max(vout*2,na.rm=TRUE)"
        } 
        if( alpha>0 ) {
            alphameth <- "sum(rep(vout^alpha,2)*(1/(window)^4),na.rm=TRUE)^(1/alpha)"
        }
        if( alpha>100 ) warning("With this alpha value you may get integer overflow: consider decreasing it.")

# Check if there are NAs in the matrices
if ( methods::is(rasterm[[1]],"SpatRaster") ){
    if( any(sapply(lapply(unlist(rasterm),length),is.na)==TRUE) )
    warning("\n One or more SpatRasters contain NAs which will be treated as 0s")
    } else if ( methods::is(rasterm[[1]],"matrix") ){
        if(any(sapply(rasterm, is.na)==TRUE) ) {
            warning("\n One or more matrices contain NAs which will be treated as 0s")
        }
    }
# Check whether the chosen distance metric is valid or not
if( dist_m=="euclidean" | dist_m=="manhattan" | dist_m=="canberra" | dist_m=="minkowski" | dist_m=="mahalanobis" ) {
## Decide what function to use
if( dist_m=="euclidean") {
    distancef <- get(".meuclidean")
    } else if( dist_m=="manhattan" ) {
        distancef <- get(".mmanhattan")
        } else if( dist_m=="canberra" ) {
            distancef <- get(".mcanberra")
            } else if( dist_m=="minkowski" ) {
                if( lambda==0 ) {
                    stop("The Minkowski distance for lambda = 0 is infinity; please choose another value for lambda.")
                    } else {
                        distancef <- get(".mminkowski") 
                    }
                    } else if( dist_m=="mahalanobis" ) {
                        distancef <- get(".mmahalanobis")
                        warning("Multimahalanobis distance is not fully supported...")
                    } 
                    } else if (is.matrix(dist_m)) {
                        distancef=dist_m
                        } else {
                            stop("Distance function not defined for multidimensional Rao's Q; please choose among euclidean, manhattan, canberra, minkowski, mahalanobis...")
                        }
# Derive Rao
tw <- apply(crop1dt, 2,function(x) {
    y <- summary(as.factor(x),maxsum=10000)
    if( "NA's"%in%names(y) ) {
        y <- y[-length(y)]
    }
    return(y)
    })

vcomb <- utils::combn(length(tw[[1]]),2)
vout <- c()
for( p in 1:ncol(vcomb) ) {
    lpair <- list(cbind(vcomb[1,p],vcomb[2,p]))
    vout[p] <- distancef(lpair)/mfactor
}
# Evaluate the parsed alpha method
mpaRaoOareaS <- eval(parse(text=alphameth))
gc()
}
return(mpaRaoOareaS)
}