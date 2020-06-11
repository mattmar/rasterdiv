mpaRaoS <- function(x,rasterm,alpha,w,dist_m,na.tolerance,rescale,lambda,diag,debugging,mfactor)
{
    mfactor <- 1 #Temporary patch
    window = 2*w+1
    diagonal <- ifelse(diag==TRUE,0,NA)
    raoqe <- matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
# Check if there are NAs in the matrices
    if ( is(x[[1]],"RasterLayer") ){
        if(any(sapply(lapply(unlist(x),length),is.na)==TRUE))
            message("\n Warning: One or more RasterLayers contain NA's which will be treated as 0")
    } else if ( is(x[[1]],"matrix") ){
        if(any(sapply(x, is.na)==TRUE) ) {
            message("\n Warning: One or more matrices contain NA's which will be treated as 0")
        }
    }
#
## Check whether the chosen distance metric is valid or not
#
    if( dist_m=="euclidean" | dist_m=="manhattan" | dist_m=="canberra" | dist_m=="minkowski" | dist_m=="mahalanobis" ) {
#
## Define distance functions
#
#euclidean
        multieuclidean <- function(x) {
            tmp <- lapply(x, function(y) {
                (y[[1]]-y[[2]])^2
            })
            return(sqrt(Reduce(`+`,tmp)))
        }
#manhattan
        multimanhattan <- function(x) {
            tmp <- lapply(x, function(y) {
                abs(y[[1]]-y[[2]])
            })
            return(Reduce(`+`,tmp))
        }
#canberra
        multicanberra <- function(x) {
            tmp <- lapply(x, function(y) {
                abs(y[[1]] - y[[2]]) / (abs(y[[1]]) + abs(y[[2]]))
            })
            return(Reduce(`+`,tmp))
        }
#minkowski
        multiminkowski <- function(x) {
            tmp <- lapply(x, function(y) {
                abs((y[[1]]-y[[2]])^lambda)
            })
            return(Reduce(`+`,tmp)^(1/lambda))
        }
#mahalanobis
        multimahalanobis <- function(x){
            tmp <- matrix(unlist(lapply(x,function(y) as.vector(y))),ncol=2)
            tmp <- tmp[!is.na(tmp[,1]),] 
            if( length(tmp)==0 | is.null(dim(tmp)) ) {
                return(NA)
            } else if(rcond(cov(tmp)) <= 0.001) {
                return(NA)
            } else {
#return the inverse of the covariance matrix of tmp; aka the precision matrix
                inverse<-solve(cov(tmp)) 
                if(debugging){
                    print(inverse)
                }
                tmp<-scale(tmp,center=T,scale=F)
                tmp<-as.numeric(t(tmp[1,])%*%inverse%*%tmp[1,])
                return(sqrt(tmp))
            }
        }
#
## Decide what function to use
#
        if( dist_m=="euclidean") {
            distancef <- get("multieuclidean")
        } else if( dist_m=="manhattan" ) {
            distancef <- get("multimanhattan")
        } else if( dist_m=="canberra" ) {
            distancef <- get("multicanberra")
        } else if( dist_m=="minkowski" ) {
            if( lambda==0 ) {
                stop("The Minkowski distance for lambda = 0 is infinity; please choose another value for lambda.")
            } else {
                distancef <- get("multiminkowski") 
            }
        } else if( dist_m=="mahalanobis" ) {
            distancef <- get("multimahalanobis")
            warning("Multimahalanobis distance is not fully supported...")
        }
    } else {
        stop("Distance function not defined for multidimensional Rao's Q; please choose among euclidean, manhattan, canberra, minkowski, mahalanobis!")
    }
    if(debugging) {
        message("#check: After distance calculation in multimenional clause.")
        print(distancef)
    }
#
## Reshape values
#
    vls<-lapply(x, function(x) {raster::as.matrix(x)})
#
## Rescale and add additional columns and rows for moving w
#
    hor<-matrix(NA,ncol=dim(vls[[1]])[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(vls[[1]])[1]+w*2)
    if(rescale) {
        trastersm<-lapply(vls, function(x) {
            t1 <- raster::scale(raster(cbind(ver,rbind(hor,x,hor),ver)))
            t2 <- raster::as.matrix(t1)
            return(t2)
        })
    } else {
        trastersm<-lapply(vls, function(x) {
            cbind(ver,rbind(hor,x,hor),ver)
        })
    }
    if(debugging) {
        message("#check: After rescaling in multimensional clause.")
        print(distancef)
    }
#
## Loop over all the pixels in the matrices
#
    if( (ncol(vls[[1]])*nrow(vls[[1]]))> 10000) {
        message("\n Warning: ",ncol(vls[[1]])*nrow(vls[[1]])*length(vls), " cells to be processed, it may take some time... \n")
    }
    for (cl in (1+w):(dim(vls[[1]])[2]+w)) {
        for(rw in (1+w):(dim(vls[[1]])[1]+w)) {
            if( length(!which(!trastersm[[1]][c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
                raoqe[rw-w,cl-w] <- NA
            } else {
                tw <- lapply(trastersm, function(x) { x[(rw-w):(rw+w),(cl-w):(cl+w)]
            })
#
## Vectorize the matrices in the list and calculate
#Among matrix pairwase distances
                lv <- lapply(tw, function(x) {as.vector(t(x))})
                vcomb <- combn(length(lv[[1]]),2)
                vout <- c()
                for(p in 1:ncol(vcomb) ) {
                    lpair <- lapply(lv, function(chi) {
                        c(chi[vcomb[1,p]],chi[vcomb[2,p]])
                    })
                    vout[p] <- distancef(lpair)
                }
                raoqe[rw-w,cl-w] <- (sum(rep(vout^alpha,2) * (1/(window)^4),na.rm=TRUE) ^ (1/alpha)) / mfactor
            }
        }
        #svMisc::progress(value=cl/(ncol(trastersm)-1)*100, max.value=100, progress.bar = F,init=T)
    }
    message(paste("\nCalculation of Multidimensional parametric Rao's index with alpha",alpha,"is completed.\n"))
    return(raoqe)
}