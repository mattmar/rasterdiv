paRaoP <- function(rasterm,alpha,w,dist_m,na.tolerance,diag,debugging,isfloat,mfactor,np) {
    # Some initial housekeeping
    message("\n\nProcessing alpha: ",alpha, " Moving Window: ", 2*w+1)
    mfactor <- ifelse(isfloat,mfactor,1) 
    window = 2*w+1
    diagonal <- ifelse(diag==TRUE,0,NA)
    # Set a progress bar
    pb <- txtProgressBar(title = "Iterative training", min = w, max = dim(rasterm)[2]+w, style = 3)
    # If alpha ~ +infinite
    if( alpha >= .Machine$integer.max | is.infinite(alpha) ) {
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
     d1 <- stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
 }
 out <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    if(debugging) {cat(paste(cl))}
    # Update progress bar
    setTxtProgressBar(pb, cl)
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
     d1<-stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
 }
 out <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    if(debugging) {
        cat(paste(cl))
    }
    # Update progress bar
    setTxtProgressBar(pb, cl)
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
                p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod,na.rm=TRUE))
                d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
                vv <- (sum(p1*(d2^alpha)*2,na.rm=TRUE)^(1/alpha)) / mfactor
                return(vv)
            }
        }
    })
    return(paRaoOP)
    Sys.sleep(5)
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
        d1<-stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
    }
    out <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
        if(debugging) {
            cat(paste(cl))
        }
        # Update progress bar
        setTxtProgressBar(pb, cl)
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
                    p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod,na.rm=TRUE))
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
