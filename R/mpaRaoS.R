mpaRaoS <- function(x,alpha,w,dist_m,na.tolerance,rescale,lambda,diag,debugging,isfloat,mfactor) {
    # Some initial housekeeping
    message("\n\nProcessing alpha: ",alpha, " Moving Window: ", 2*w+1)
    mfactor <- ifelse(isfloat,mfactor,1) 
    window = 2*w+1
    diagonal <- ifelse(diag==TRUE,0,NA)
    rasterm <- x[[1]]
    # Evaluate Rao's method given alpha
    if( alpha >= .Machine$integer.max | is.infinite(alpha) ) {
        alphameth <- "max(vout*2,na.rm=TRUE)"
    } else if( alpha>0 ) {
        if( alpha >100 ) warning("With this alpha value you may get integer overflow. Consider decreasing the value of alpha.")
        alphameth <- "sum(rep(vout^alpha,2)*(1/(window)^4),na.rm=TRUE)^(1/alpha)"
    } else if( alpha==0 ) {
        alphameth <- "prod(vout,na.rm=TRUE)^(1/(window^4))"
    }
    # Set a progress bar
    pb <- progress_bar$new(
        format = "\n [:bar] :elapsed -- Approximate ETA: :eta \n",
        total = (dim(rasterm)[2]+w), 
        clear = FALSE, 
        width = 80, 
        force = TRUE)
    # Define output matrix
    raoqe <- matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
    # Check if there are NAs in the matrices
    if ( is(x[[1]],"RasterLayer") ){
        if(any(sapply(lapply(unlist(x),length),is.na)==TRUE))
            warning("\n One or more RasterLayers contain NAs which will be treated as 0s")
    } else if ( is(x[[1]],"matrix") ){
        if(any(sapply(x, is.na)==TRUE) ) {
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
    } else {
        stop("Distance function not defined for multidimensional Rao's Q; please choose among euclidean, manhattan, canberra, minkowski, mahalanobis...")
    }
    if(debugging) {
        message("#check: After distance calculation in multidimensional function.")
    }
    # Add additional columns and rows to account for moving window size
    hor <- matrix(NA,ncol=dim(x[[1]])[2],nrow=w)
    ver <- matrix(NA,ncol=w,nrow=dim(x[[1]])[1]+w*2)
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
    for (cl in (1+w):(dim(x[[1]])[2]+w)) {
        # Update progress bar
        pb$tick()
        for(rw in (1+w):(dim(x[[1]])[1]+w)) {
            if( length(!which(!trastersm[[1]][c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) <= (window^2-((window^2)*na.tolerance)) ) {
                raoqe[rw-w,cl-w] <- NA
            } else {
                tw <- lapply(trastersm, function(x) { 
                    x[(rw-w):(rw+w),(cl-w):(cl+w)]
                })
                ## Vectorise the matrices in the list and calculate among matrix pairwase distances
                lv <- lapply(tw, function(x) {as.vector(t(x))})
                vcomb <- combn(length(lv[[1]]),2)
                vout <- c()
                # Exclude windows with only 1 category
                if( sum(sapply(lv, function(x) length(unique(x))))<3 ) {
                    raoqe[rw-w,cl-w] <- 0
                } else {
                    for( p in 1:ncol(vcomb) ) {
                        lpair <- lapply(lv, function(chi) {
                            c(chi[vcomb[1,p]],chi[vcomb[2,p]])
                        })
                        vout[p] <- distancef(lpair)/mfactor
                    }
                }
                # Evaluate the parsed alpha method
                raoqe[rw-w,cl-w] <- eval(parse(text=alphameth))
            }
        }
    }
    return(raoqe)
}