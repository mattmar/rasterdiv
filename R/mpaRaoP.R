#' Multidimensional parallel Parametric Rao's index of quadratic entropy (Q)
#'
#' Multidimensional parametric Rao's index of quadratic entropy (Q).
#'
#' @param x input list.
#' @param alpha alpha value for order of diversity in Hill's Index.
#' @param window half of the side of the square moving window.
#' @param dist_m Type of distance used.
#' @param na.tolerance a numeric value \eqn{(0.0-1.0)} which indicates the proportion
#'   of NA values that will be tolerated to calculate Rao's index in each moving
#'   window over \emph{x}. If the relative proportion of NA's in a moving window is
#'   bigger than na.tolerance, then the value of the window will be set as NA,
#'   otherwise Rao's index will be calculated considering the non-NA values.
#'   Default values is 0.0 (i.e., no tolerance for NA's).
#' @param rescale Scale and centre values in each of the element of x.
#' @param lambda Lambda value for Minkowski distance.
#' @param diag Boolean. Diagonal of the distance matrix.
#' @param time_vector time; 
#' @param stepness numeric; steepness of the logistic function.
#' @param midpoint numeric; midpoint of the logistic function
#' @param cycle_length string; The length of the cycle. Can be a numeric value or a string specifying the units ('year', 'month', 'day', 'hour', 'minute', 'second'). When numeric, the cycle length is in the same units as time_scale. When a string, it specifies the time unit of the cycle.
#' @param time_scale string; Specifies the time scale for the conversion. Must be one of 'year', 'month', 'day', 'hour', 'minute', 'second'. When cycle_length is a string, time_scale changes the unit in which the result is expressed. When cycle_length is numeric, time_scale is used to compute the elapsed time in seconds.
#' @param debugging a boolean variable set to FALSE by default. If TRUE, additional
#'   messages will be printed. For de-bugging only.
#' @param isfloat Are the input data floats?
#' @param mfactor Multiplication factor in case of input data as float numbers.
#' @param np the number of processes (cores) which will be spawned.
#' @param progBar logical. If TRUE a progress bar is shown.
#'
#' @return A list of matrices of dimension \code{dim(x)} with length equal to the
#'   length of \code{alpha}.
#'
#' @author Duccio Rocchini \email{duccio.rocchini@unibo.it}, Marcantonio Matteo
#'   \email{marcantoniomatteo@gmail.com}
#'
#' @seealso \code{\link{paRao}}
#'
#' @keywords internal

mpaRaoP <- function(x,alpha,window,dist_m,na.tolerance,rescale,lambda, diag, time_vector, stepness, midpoint, cycle_length, time_scale, debugging, isfloat, mfactor, np, progBar) {
   # `win` is the operative moving window
   win = window 
   NAwin <- 2*window+1
   message("\n\nProcessing alpha: ",alpha, " Moving Window: ", NAwin)
    
    # Set a progress bar
    if( progBar ) {
        pb <- progress::progress_bar$new(
        format = "[:bar] :percent in :elapsed\n",
        # Total number of ticks is the number of column +NA columns divided the number of processor.
        total = (dim(x[[1]])[2]/np)+5, 
        clear = FALSE, 
        width = 60, 
        force = FALSE)
    }

    mfactor <- ifelse(isfloat,mfactor,1) 
    diagonal <- ifelse(diag==TRUE,0,NA)
    rasterm <- x[[1]]
    # Evaluate Rao's method given alpha
    if( (alpha>=.Machine$integer.max) | is.infinite(alpha) ) {
        alphameth <- "max(vout*2,na.rm=TRUE)"
        } else if( alpha>0 ) {
            if( alpha >100 ) warning("With this alpha value you may get integer overflow. Consider decreasing the value of alpha.")
            alphameth <- "sum((rep(vout^alpha,2) * (1/(NAwin)^4)),na.rm=TRUE) ^ (1/alpha)"
            } else if( alpha==0 ) {
                alphameth <- "prod(vout,na.rm=TRUE) ^ (1/(NAwin^4))"
                } else {
                    stop()
                }
    # Check if there are NAs in the matrices
    if ( methods::is(x[[1]],"SpatRaster") ){
        if(any(sapply(lapply(unlist(x),length),is.na)==TRUE))
        warning("\n One or more SpatRasters contain NA's which will be treated as 0")
        } else if ( methods::is(x[[1]],"matrix") ){
            if(any(sapply(x, is.na)==TRUE) ) {
                warning("\n One or more matrices contain NA's which will be treated as 0")
            }
        }
# Validate and set the distance function
validDistanceMetrics <- c("euclidean", "manhattan", "canberra", "minkowski", "mahalanobis", "twdtw")
if (dist_m %in% validDistanceMetrics) {
    switch(dist_m,
        euclidean = distancef <- get(".meuclidean"),
        manhattan = distancef <- get(".mmanhattan"),
        canberra = distancef <- get(".mcanberra"),
        twdtw = distancef <- get(".mtwdtw"),
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
            stop("Invalid distance metric. Choose among 'euclidean', 'manhattan', 'canberra', 'minkowski', 'mahalanobis', 'twdtw' or provide a matrix.")
        }
        # Debugging check
        if (debugging) {
            message("#check: After setting up distance calculation in multidimensional Rao's Q function.")
        }
    # Add additional columns and rows to account for moving NAwin size
    hor <- matrix(NA,ncol=dim(x[[1]])[2],nrow=win)
    ver <- matrix(NA,ncol=win,nrow=dim(x[[1]])[1]+win*2)
    trastersm <- lapply(x, function(x) {
        cbind(ver,rbind(hor,x,hor),ver)
        })
    if(debugging) {
        message("#check: After rescaling in multimensional clause.")
        print(distancef)
    }
    # Loop over all the pixels in the matrices
    if( (ncol(x[[1]])*nrow(x[[1]]))>10000 ) {
        warning("",ncol(x[[1]])*nrow(x[[1]])*length(x), " cells process, it may take quite some time... \n")
    }
    # Parallelised parametric multidimensional Rao
    out <- foreach::foreach(cl=(1+win):(dim(rasterm)[2]+win),.verbose = F, .export=c("alpha")) %dopar% {
        # Update progress bar
        if(progBar) pb$tick()
        # Row loop
        mpaRaoOP <- sapply((1+win):(dim(rasterm)[1]+win), function(rw) {
            if(debugging) {
                message("#check: Inside sapply.")
            }
            if( length(!which(!trastersm[[1]][c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) < floor(NAwin^2-((NAwin^2)*na.tolerance)) ) {
                vv <- NA
                return(vv)
                } else {
                    tw <- lapply(trastersm, function(x) { 
                        x[(rw-win):(rw+win),(cl-win):(cl+win)]
                        })
                # Vectorise the matrices in the list and calculate between matrices pairwase distances
                lv <- lapply(tw, function(x) as.vector(t(x)))
                vcomb <- utils::combn(length(lv[[1]]),2)
                # Exclude windows with only 1 category in all lists
                if( sum(sapply(lv, function(x) length(unique(x))),na.rm=TRUE)<(length(lv)+1) ) {
                    vv <- 0
                    } else {
                        vout <- sapply(1:ncol(vcomb), function(p) {
                            lpair <- lapply(lv, function(chi) {
                                c(chi[vcomb[1,p]],chi[vcomb[2,p]])
                                })
                            return(
                                if (dist_m == "twdtw") {
                                    llist <- list(sapply(lpair, function(x) x[1]), sapply(lpair, function(x) x[2]))
                                    distancef(llist, time_vector = time_vector, stepness = stepness, midpoint = midpoint, cycle_length = cycle_length, time_scale = time_scale) / mfactor
                                    } else {
                                        distancef(lpair) / mfactor
                                        })
                            })
                            # Evaluate the parsed alpha method
                            vv <- eval(parse(text=alphameth))
                        }
                        return(vv)
                    }
                    })
        return(mpaRaoOP)
    }
    return(do.call(cbind,out))
}