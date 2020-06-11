paRaoS <- function(rasterm,alpha,w,dist_m,na.tolerance,diag,debugging,isfloat,mfactor) 
{
	message("\n\nProcessing alpha ",alpha)
	mfactor <- ifelse(isfloat,mfactor,1) 
	window = 2*w+1
	diagonal <- ifelse(diag==TRUE,0,NA)
	# If alpha ~ +infinite
	if( alpha >= .Machine$integer.max ) {
		paRaoOS <- matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
		# Reshape values
		values <- as.numeric(as.factor(rasterm))
		rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
		# Add additional columns and rows for moving window
		hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
		ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
		trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
		# Derive distance matrix
		classes <- levels(as.factor(rasterm))
		if( is.character(dist_m) | is.function(dist_m) ) {
			d1<-proxy::dist(as.numeric(classes),method=dist_m)
		} else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
			d1<-stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
		}
		# Loop over each pixel
		for (cl in (1+w):(dim(rasterm)[2]+w)) {
			for(rw in (1+w):(dim(rasterm)[1]+w)) {
				if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
					paRaoOS[rw-w,cl-w]<-NA
				} else {
					tw<-summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
					if( "NA's"%in%names(tw) ) {
						tw <- tw[-length(tw)]
					}
					if(debugging) {
						message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window^2)
					}
					tw_labels <- names(tw)
					tw_values <- as.vector(tw)
          			# if clause to exclude windows with only 1 category
					if(length(tw_values) < 2) {
						paRaoOS[rw-w,cl-w] <- NA
					} else {
						p <- tw_values/sum(tw_values)
						p1 <- diag(diagonal,length(tw_values))
						p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod,na.rm=TRUE))
						d2 <- unname(proxy::as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
						paRaoOS[rw-w,cl-w] <- max(d2*2,na.rm=TRUE) / mfactor
					}
				}
			} 
		}
		return(paRaoOS)
	# If alpha > 0
	} else if( alpha>0 ){
		paRaoOS<-matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
		# Reshape values
		values<-as.numeric(as.factor(rasterm))
		rasterm_1<-matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
		# Add additional columns and rows for moving window
		hor<-matrix(NA,ncol=dim(rasterm)[2],nrow=w)
		ver<-matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
		trasterm<-cbind(ver,rbind(hor,rasterm_1,hor),ver)
		# Derive distance matrix
		classes<-levels(as.factor(rasterm))
		if( is.character(dist_m) | is.function(dist_m) ) {
			d1<-proxy::dist(as.numeric(classes),method=dist_m)
		} else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
			d1<-stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
		}
		# Loop over each pixel
		for(cl in (1+w):(dim(rasterm)[2]+w)) {
			for(rw in (1+w):(dim(rasterm)[1]+w)) {
				if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
					paRaoOS[rw-w,cl-w]<-NA
				} else {
					tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
					if( "NA's"%in%names(tw) ) {
						tw <- tw[-length(tw)]
					}
					if(debugging) {
						message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window^2)
					}
					tw_labels <- names(tw)
					tw_values <- as.vector(tw)
        			# if clause to exclude windows with only 1 category
					if( length(tw_values) < 2 ) {
						paRaoOS[rw-w,cl-w] <- NA
					} else {
						p <- tw_values/sum(tw_values,na.rm=TRUE)
						p1 <- diag(diagonal,length(tw_values))
						p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod,na.rm=TRUE))
						d2 <- unname(proxy::as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
						paRaoOS[rw-w,cl-w] <- ( (sum((p1)*(d2^alpha)*2,na.rm=TRUE))^(1/alpha) ) / mfactor
					}
				}			} 
			}
			#message(("\nCalculation of sequential Parametric Rao's index complete.\n"))
			return(paRaoOS)
		# If alpha == 0
		} else if( alpha==0 ){
			paRaoOS <- matrix(rep(NA,dim(rasterm)[1]*dim(rasterm)[2]),nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
			# Reshape values
			values <- as.numeric(as.factor(rasterm))
			rasterm_1 <- matrix(data=values,nrow=dim(rasterm)[1],ncol=dim(rasterm)[2])
			# Add additional columns and rows for moving window
			hor <- matrix(NA,ncol=dim(rasterm)[2],nrow=w)
			ver <- matrix(NA,ncol=w,nrow=dim(rasterm)[1]+w*2)
			trasterm <- cbind(ver,rbind(hor,rasterm_1,hor),ver)
			# Derive distance matrix
			classes <- levels(as.factor(rasterm))
			if( is.character(dist_m) | is.function(dist_m) ) {
				d1 <- proxy::dist(as.numeric(classes),method=dist_m)
			} else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
				d1 <- stats::as.dist(xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
			}
			for (cl in (1+w):(dim(rasterm)[2]+w)) {
				for(rw in (1+w):(dim(rasterm)[1]+w)) {
					if( length(!which(!trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]%in%NA)) < window^2-((window^2)*na.tolerance) ) {
						paRaoOS[rw-w,cl-w] <- NA
					} else {
						tw <- summary(as.factor(trasterm[c(rw-w):c(rw+w),c(cl-w):c(cl+w)]),maxsum=10000)
						if( "NA's"%in%names(tw) ) {
							tw <- tw[-length(tw)]
						}
						if(debugging) {
							message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window^2)
						}
						tw_labels <- names(tw)
						tw_values <- as.vector(tw)
      				# if clause to exclude windows with only 1 category
						if(length(tw_values) < 2) {
							paRaoOS[rw-w,cl-w] <- NA
						} else {
							d2 <- unname( proxy::as.matrix(d1,diag=diagonal)[as.numeric(tw_labels),as.numeric(tw_labels)] )
							paRaoOS[rw-w,cl-w] <- ( prod(d2/mfactor,na.rm=TRUE)^(1/(window^2)) )
						}
					} 
				} 
			}
			return(paRaoOS)
			#message(("\nCalculation of sequential Parametric Rao's index complete.\n"))
		}
	}