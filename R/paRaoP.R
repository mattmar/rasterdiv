paRaoP <- function(rasterm,alpha,w,dist_m,na.tolerance,diag,debugging,isfloat,mfactor) 
{
    #Some initial tricks
    message("\n\nProcessing alpha ",alpha)
    mfactor <- ifelse(isfloat,mfactor,1) 
    window = 2*w+1
    diagonal <- ifelse(diag==TRUE,0,NA)
    # If alpha ~ +infinite
    if( alpha >= .Machine$integer.max ) {
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
   } else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
       d1<-stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
   }
   out <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    if(debugging) {
        cat(paste(cl))
    }
    paRaoOP <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
        if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
            vv<-NA
            return(vv)
        } 
        else {
            tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            if( "NA's"%in%names(tw) ) {
                tw<-tw[-length(tw)]
            }
            if( debugging ) {
                message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window^2)
            }
            tw_labels <- names(tw)
            tw_values <- as.vector(tw)
                #if clause to exclude windows with only 1 category
            if( length(tw_values) <2 ) {
                vv<-NA
                return(vv)
            }
            else {
                p <- tw_values/sum(tw_values)
                p1 <- diag(0,length(tw_values))
                p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod,na.rm=TRUE))
                d2 <- unname(as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
                vv <- max(d2*2,na.rm=TRUE) / mfactor
                return(vv)
            }
        }
    })
    return(paRaoOP)
} #End classic Parametric Rao - parallelized
    #message(("\n\nCalculation of Parametric Rao's index complete.\n"))
return(do.call(cbind,out))
} else if( alpha>0 ) {
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
   } else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
       d1<-stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
   }
   out <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    if(debugging) {
        cat(paste(cl))
    }
    paRaoOP <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
        if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
            vv<-NA
            return(vv)
        } 
        else {
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
                vv<-NA
                return(vv)
            }
            else {
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
} #End classic Parametric Rao - parallelized
return(do.call(cbind,out))
} else if( alpha==0 ) {
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
   } else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
       d1<-stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
   }
   out <- foreach(cl=(1+w):(dim(rasterm)[2]+w),.verbose = F) %dopar% {
    if(debugging) {
        cat(paste(cl))
    }
    paRaoOP <- sapply((1+w):(dim(rasterm)[1]+w), function(rw) {
        if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
            vv<-NA
            return(vv)
        } 
        else {
            tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
            if( "NA's"%in%names(tw) ) {
                tw<-tw[-length(tw)]
            }
            if( debugging ) {
                message("Working on coords ",rw,",",cl,". classes length: ",length(tw),". window size=",window^2)
            }
            tw_labels <- names(tw)
            tw_values <- as.vector(tw)
            #if clause to exclude windows with only 1 category
            if( length(tw_values) <2 ) {
                vv<-NA
                return(vv)
            }
            else {
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
} #End classic Parametric Rao - parallelized
return(do.call(cbind,out))
}
}
