mpaRaoAreaS <- function(rasterm, area, alpha, simplify, dist_m, rescale, lambda, window) {
# Some initial housekeeping
if(alpha<=0) {
    stop("Alpha<=0 not yet implemented for areas.")
}
mfactor <- 100^simplify
crop1 <- crop(rasterm, area)
crop1dt <- raster::as.matrix(crop1)*mfactor
storage.mode(crop1dt) <- "integer"
classes <- levels(as.factor(crop1dt))

# Evaluate Rao's method given alpha
window <- nrow(crop1dt) #Temporary patch?
if( alpha >= .Machine$integer.max | is.infinite(alpha) ) {
    alphameth <- "max(vout*2,na.rm=TRUE)"
} 
if( alpha>0 ) {
    alphameth <- "sum(rep(vout^alpha,2)*(1/(window)^4),na.rm=TRUE)^(1/alpha)"
}
if( alpha>100 ) warning("With this alpha value you may get integer overflow. Consider decreasing the value of alpha.")

# Check if there are NAs in the matrices
if ( is(rasterm[[1]],"RasterLayer") ){
    if(any(sapply(lapply(unlist(rasterm),length),is.na)==TRUE))
    warning("\n One or more RasterLayers contain NAs which will be treated as 0s")
    } else if ( is(rasterm[[1]],"matrirasterm") ){
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

vcomb <- combn(nrow(tw),2)
vout <- c()
for( p in 1:ncol(vcomb) ) {
    lpair <- list(cbind(vcomb[1,p],vcomb[2,p]))
    vout[p] <- distancef(lpair)/mfactor
}
# Evaluate the parsed alpha method
mpaRaoOareaS <- eval(parse(text=alphameth))
gc()
return(mpaRaoOareaS)
}