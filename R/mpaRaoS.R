#' Multidimensional sequential Parametric Rao's index of quadratic entropy (Q)
#'
#' This function calculates the multidimensional parametric Rao's index of quadratic entropy (Q) using a 
#' sequential method. It is particularly useful in contexts where parallel computation is not feasible or desired.
#' The function applies a moving window approach to the provided raster data stack.
#'
#' @param x input list.
#' @param alpha Numeric; alpha value for order of diversity in Hill's Index.
#' @param window Numeric; half of the side of the square moving window used for calculation.
#' @param dist_m Character; type of distance used in the analysis.
#' @param na.tolerance Numeric; a threshold between 0.0 and 1.0 indicating the allowable proportion of NA 
#' values within each moving window. If the proportion of NA values exceeds this, the window's value is set as 
#' NA; otherwise, the computation uses the non-NA values.
#' @param rescale Logical; if TRUE, scales and centres the values in each element of 'x'.
#' @param lambda Numeric; lambda value used for Minkowski distance calculation.
#' @param diag Logical; if TRUE, includes the diagonal of the distance matrix in computations.
#' @param debugging Logical; if TRUE, additional diagnostic messages are output, useful for debugging. Default 
#' is FALSE.
#' @param isfloat Logical; specifies if the input data are floats.
#' @param mfactor Numeric; multiplication factor applied if input data are float numbers.
#' @param np Number of processes for parallel computation.
#' @return A list of matrices, each representing a layer of the input RasterStack, containing calculated 
#' Rao's index values. The dimensions correspond to those of the input, and the list length is equal to the 
#' length of 'alpha'.
#' @seealso \code{\link{paRao}} for the parallelized version of the Rao's index computation.
#' @author Duccio Rocchini \email{duccio.rocchini@@unibo.it}, 
#' Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}

mpaRaoS <- function(x, alpha, window, dist_m, na.tolerance, rescale, lambda, diag, debugging, isfloat, mfactor, np) {
    # `win` is the operative moving window
    win = window 
    NAwin <- 2*window+1
    message("\n\nProcessing alpha: ",alpha, " Moving Window: ", NAwin)
    # Set a progress bar
    pb <- progress::progress_bar$new(
        format = "[:bar] :percent in :elapsed\n",
    # Total number of ticks is the number of column +NA columns divided the number of processor.
    total = dim(x[[1]])[2], 
    clear = FALSE, 
    width = 60, 
    force = FALSE)
    
    mfactor <- ifelse(isfloat,mfactor,1) 
    diagonal <- ifelse(diag==TRUE,0,NA)
    rasterm <- x[[1]]
    # Evaluate Rao's method given alpha
    if( alpha >= .Machine$integer.max | is.infinite(alpha) ) {
        alphameth <- "max(vout*2,na.rm=TRUE)"
        } else if( alpha>0 ) {
            if( alpha >100 ) warning("With this alpha value you may get integer overflow. Consider decreasing the value of alpha.")
            alphameth <- "sum(rep(vout^alpha,2)*(1/(NAwin)^4),na.rm=TRUE)^(1/alpha)"
            } else if( alpha==0 ) {
                alphameth <- "prod(vout,na.rm=TRUE)^(1/(NAwin^4))"
            }
    # Define output matrix
    raoqe <- matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
# Check for NAs in SpatRasters or matrices
if (methods::is(x[[1]], "SpatRaster")) {
    if (any(sapply(x, function(rast) any(is.na(terra::values(rast)))))) {
        warning("One or more SpatRasters contain NAs, which will be treated as 0s.")
    }
    } else if (methods::is(x[[1]], "matrix")) {
        if (any(sapply(x, is.na))) {
            warning("One or more matrices contain NAs, which will be treated as 0s.")
        }
    }

# Validate and set the distance function
validDistanceMetrics <- c("euclidean", "manhattan", "canberra", "minkowski", "mahalanobis")
if (dist_m %in% validDistanceMetrics) {
    switch(dist_m,
        euclidean = distancef <- get(".meuclidean"),
        manhattan = distancef <- get(".mmanhattan"),
        canberra = distancef <- get(".mcanberra"),
        minkowski = {
            if (lambda == 0) stop("Minkowski distance with lambda = 0 is undefined. Choose another value.")
            distancef <- get(".mminkowski")
            },
            mahalanobis = {
                distancef <- get(".mmahalanobis")
                warning("Mahalanobis distance is not fully supported for multidimensional Rao's Q.")
            }
            )
    } else if (is.matrix(dist_m)) {
        distancef <- dist_m
        } else {
            stop("Invalid distance metric. Choose among 'euclidean', 'manhattan', 'canberra', 'minkowski', 'mahalanobis', or provide a matrix.")
        }
        # Debugging check
        if (debugging) {
            message("#check: After setting up distance calculation in multidimensional Rao's Q function.")
        }
    # Add additional columns and rows to account for moving window size
    hor <- matrix(NA,ncol=dim(x[[1]])[2],nrow=win)
    ver <- matrix(NA,ncol=win,nrow=dim(x[[1]])[1]+win*2)
    trastersm <- lapply(x, function(x) {
        cbind(ver,rbind(hor,x,hor),ver)
        })
    if(debugging) {
        message("#check: After adding columns in multimensional function.")
        print(distancef)
    }
    ## Loop over all the pixels in the matrices
    if( (ncol(x[[1]])*nrow(x[[1]]))>10000 ) {
        message("\n Warning: ",ncol(x[[1]])*nrow(x[[1]])*length(x), " cells to be processed, it may take some time... \n")
    }
    for (cl in (1+win):(dim(x[[1]])[2]+win)) {
        # Update progress bar
        pb$tick()
        for(rw in (1+win):(dim(x[[1]])[1]+win)) {
            if( length(!which(!trastersm[[1]][c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
                raoqe[rw-win,cl-win] <- NA
                } else {
                    tw <- lapply(trastersm, function(x) { 
                        x[(rw-win):(rw+win),(cl-win):(cl+win)]
                        })
                ## Vectorise the matrices in the list and calculate among matrix pairwase distances
                lv <- lapply(tw, function(x) {as.vector(t(x))})
                vcomb <- utils::combn(length(lv[[1]]),2)
                vout <- c()
                # Exclude windows with only 1 category
                if( sum(sapply(lv, function(x) length(unique(x))))<3 ) {
                    raoqe[rw-win,cl-win] <- 0
                    } else {
                        for( p in 1:ncol(vcomb) ) {
                            lpair <- lapply(lv, function(chi) {
                                c(chi[vcomb[1,p]],chi[vcomb[2,p]])
                                })
                            vout[p] <- distancef(lpair)/mfactor
                        }
                    }
                # Evaluate the parsed alpha method
                raoqe[rw-win,cl-win] <- eval(parse(text=alphameth))
            }
        }
    }
    return(raoqe)
}