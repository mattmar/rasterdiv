# Distance function for multidimension Rao
# euclidean
.meuclidean <- function(x) {
    tmp <- lapply(x, function(y) {
        (y[[1]]-y[[2]])^2
    })
    return(sqrt(Reduce(`+`,tmp)))
}
# manhattan
.mmanhattan <- function(x) {
    tmp <- lapply(x, function(y) {
        abs(y[[1]]-y[[2]])
    })
    return(Reduce(`+`,tmp))
}
# canberra
.mcanberra <- function(x) {
    tmp <- lapply(x, function(y) {
        abs(y[[1]] - y[[2]]) / (abs(y[[1]]) + abs(y[[2]]))
    })
    return(Reduce(`+`,tmp))
}
# minkowski
.mminkowski <- function(x, lambda=lambda) {
    tmp <- lapply(x, function(y) {
        abs((y[[1]]-y[[2]])^lambda)
    })
    return(Reduce(`+`,tmp)^(1/lambda))
}
# mahalanobis
.mmahalanobis <- function(x, debugging=debugging){
    tmp <- matrix(unlist(lapply(x,function(y) as.vector(y))),ncol=2)
    tmp <- tmp[!is.na(tmp[,1]),] 
    if( length(tmp)==0 | is.null(dim(tmp)) ) {
        return(NA)
    } else if(rcond(stats::cov(tmp)) <= 0.001) {
        return(NA)
    } else {
        # return the inverse of the covariance matrix of tmp; aka the precision matrix
        inverse <- solve(stats::cov(tmp)) 
        if(debugging){
            print(inverse)
        }
        tmp<-scale(tmp,center=T,scale=F)
        tmp<-as.numeric(t(tmp[1,])%*%inverse%*%tmp[1,])
        return(sqrt(tmp))
    }
}