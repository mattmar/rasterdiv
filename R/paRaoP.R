#' Parallelized Parametric Rao's index of quadratic entropy (Q)
#'
#' This function computes the parametric Rao's index of quadratic entropy (Q), a measure of biodiversity 
#' that considers the evolutionary distances between species, utilizing parallel computing for enhanced 
#' performance. The computation is applied over a moving window across the input data.
#'
#' @param rasterm Matrix or data frame; the input data over which the index calculation is performed.
#' @param alpha Numeric; specifies the alpha value for the order of diversity in Hill's Index.
#' @param w Numeric; half of the side length of the square moving window used in the calculation.
#' @param dist_m Character; specifies the type of distance metric used in calculations.
#' @param na.tolerance Numeric; the threshold proportion of NA values allowed in the moving window. 
#' If exceeded, the calculation for that window is skipped. Values range from 0.0 (no tolerance) to 1.0.
#' @param diag Logical; indicates whether the diagonal of the distance matrix should be included in the 
#' computation. Typically set to FALSE.
#' @param debugging Logical; set to FALSE by default. If TRUE, additional console messages will be 
#' displayed for debugging purposes.
#' @param isfloat Logical; indicates whether the input data values are floating-point numbers.
#' @param mfactor Integer; multiplication factor in case of input data as float numbers.
#' @return A list of matrices corresponding to the computed Rao's index values. Each matrix in the list 
#' represents the calculations performed over the moving window, with dimensions equal to \code{dim(rasterm)}.
#' @author Duccio Rocchini \email{duccio.rocchini@@unibo.it},
#' Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}
#' @seealso \code{\link{paRao}} for the related non-parallelized function.

paRaoP <- function(rasterm,alpha,w,dist_m,na.tolerance,diag,debugging,isfloat,mfactor) {
    # Some initial housekeeping
    message("\n\nProcessing alpha: ",alpha, " Moving Window: ", 2*w+1)
    mfactor <- ifelse(isfloat,mfactor,1) 
    window = 2*w+1
    diagonal <- ifelse(diag==TRUE,0,NA)
    tdist <- proxy::dist(as.numeric(levels(as.factor(rasterm))),method=dist_m)
    # Min and max dist for initial checks on possible infinite or 0 operations
    maxd <- max(proxy::dist(as.numeric(levels(as.factor(rasterm))),method=dist_m))
    mind <- min(tdist[tdist>0])
    # Set up a progress bar
    pb <- utils::txtProgressBar(title = "Iterative training", min = w, max = dim(rasterm)[2]+w, style = 3)
    # If alpha ~ +infinite
    if( alpha >= .Machine$integer.max | is.infinite(alpha) | is.infinite(maxd^alpha) | (dist_m=="canberra" & mind^alpha==0) ) {
      # Reshape values
      values <- as.numeric(as.factor(rasterm))
      rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
      # Add additional columns and rows to match moving window
      hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
      ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
      trasterm <- cbind(ver,rbind(hor,rasterm_1,hor),ver)
      rm(hor,ver,rasterm_1,values)
      gc()
      if( debugging ) {
        cat("#check: Parametric Rao parallel function.")
    }
    # Derive distance matrix
    if( is.character( dist_m) | is.function(dist_m) ) {
        d1 <- proxy::dist(as.numeric(levels(as.factor(rasterm))),method=dist_m) 
    } else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
       d1 <- stats::as.dist(stats::xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
   }
   out <- foreach::foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    if(debugging) {cat(paste(cl))}
    # Update progress bar
    utils::setTxtProgressBar(pb, cl)
    # Row loop
    paRaoOP <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
        if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) <= (window^2-((window^2)*na.tolerance)) ) {
            vv <- NA
            return(vv)
        }else{
            tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            if( "NA's"%in%names(tw) ) {
                tw <- tw[-length(tw)]
            }
            if( debugging ) {
                message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window^2)
            }
            tw_labels <- names(tw)
            tw_values <- as.vector(tw)
            #if clause to exclude windows with only 1 category
            if( length(tw_values) < 2 ) {
                vv <- 0
                return(vv)
            }else{
                d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
                vv <- max(d2*2,na.rm=TRUE) / mfactor
                return(vv)
            }
        }
    })
    return(paRaoOP)
} #End classic Parametric Rao - parallelized
return(do.call(cbind,out))
# If alpha is > 0
}else if( alpha>0 ) {
    #
    ##Reshape values
    #
    values <- as.numeric(as.factor(rasterm))
    rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
    #
    ##Add additional columns and rows to match moving window
    #
    hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
    trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
    rm(hor,ver,rasterm_1,values); gc()
    if(debugging){cat("#check: Parametric Rao parallel function.")}
    #       
    ##Derive distance matrix
    #
    if( is.character( dist_m) | is.function(dist_m) ) {
       d1<-proxy::dist(as.numeric(levels(as.factor(rasterm))),method=dist_m)
   }else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
       d1<-stats::as.dist(stats::xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
   }
   out <- foreach::foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    if(debugging) {
        cat(paste(cl))
    }
    # Update progress bar
    utils::setTxtProgressBar(pb, cl)
    # Row loop
    paRaoOP <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
        if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) <= (window^2-((window^2)*na.tolerance)) ) {
            vv <- NA
            return(vv)
        }else{
            tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            if( "NA's"%in%names(tw) ) {
                tw<-tw[-length(tw)]
            }
            if( debugging ) {
                message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window^2)
            }
            tw_labels <- names(tw)
            tw_values <- as.vector(tw)
            #if clause to exclude windows with only 1 category
            if( length(tw_values) < 2 ) {
                vv <- 0
                return(vv)
            }else{
                p <- tw_values/sum(tw_values,na.rm=TRUE)
                p1 <- diag(0,length(tw_values))
                p1[lower.tri(p1)] <- c(utils::combn(p,m=2,FUN=prod,na.rm=TRUE))
                d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
                vv <- (sum((p1)*(d2^alpha)*2,na.rm=TRUE)^(1/alpha) ) / mfactor
                return(vv)
            }
        }
    })
    return(paRaoOP)
} #End classic Parametric Rao - parallelized
return(do.call(cbind,out))
}else if( alpha==0 ) {
    #
    ##Reshape values
    #
    values<-as.numeric(as.factor(rasterm))
    rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
    #
    ##Add additional columns and rows to match moving window
    #
    hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
    ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
    trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
    rm(hor,ver,rasterm_1,values); gc()
    if(debugging){cat("#check: Parametric Rao parallel function.")}
    #       
    ##Derive distance matrix
    #
    if( is.character( dist_m) | is.function(dist_m) ) {
        d1<-proxy::dist(as.numeric(levels(as.factor(rasterm))),method=dist_m)
    }else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
        d1<-stats::as.dist(stats::xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
    }
    out <- foreach::foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
        if(debugging) {
            cat(paste(cl))
        }
        # Update progress bar
        utils::setTxtProgressBar(pb, cl)
        # Row loop
        paRaoOP <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
            if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) <= (window^2-((window^2)*na.tolerance)) ) {
                vv <- NA
                return(vv)
            }else{
                tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
                if( "NA's"%in%names(tw) ) {
                    tw<-tw[-length(tw)]
                }
                if( debugging ) {
                    message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window^2)
                }
                tw_labels <- names(tw)
                tw_values <- as.vector(tw)
            #if clause to exclude windows with only 1 category
                if( length(tw_values) < 2 ) {
                    vv <- 0
                    return(vv)
                }else{
                    p <- tw_values/sum(tw_values)
                    p1 <- diag(0,length(tw_values))
                    p1[lower.tri(p1)] <- c(utils::combn(p,m=2,FUN=prod,na.rm=TRUE))
                    d2 <- unname( proxy::as.matrix(d1,diag=diagonal)[as.numeric(tw_labels),as.numeric(tw_labels)] )
                    vv <- (prod(d2,na.rm=TRUE)^(1/(window^2))) / mfactor
                    return(vv)
                }
            }
        })
        return(paRaoOP)
    }
    return(do.call(cbind,out))
}else{stop("Something went wrong. Exiting...")}
}
