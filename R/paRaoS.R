#' Sequential Parametric Rao's index of quadratic entropy (Q)
#'
#' Computes the sequential version of the parametric Rao's index of quadratic entropy (Q), 
#' a measure used in environmental and ecological studies to assess biodiversity by considering 
#' the evolutionary distance between species. The function performs calculations in a sequential 
#' manner over a moving window across the input data.
#'
#' @param x Matrix or data frame; the input data over which the index calculation is performed.
#' @param alpha Numeric; specifies the alpha value for the order of diversity in Hill's Index.
#' @param window Numeric; half of the side length of the square moving window used in the calculation.
#' @param dist_m Character; specifies the type of distance metric used in calculations.
#' @param na.tolerance Numeric; the threshold proportion of NA values allowed in the moving window. 
#' If exceeded, the calculation for that window is skipped. Values range from 0.0 (no tolerance) to 1.0.
#' @param diag Logical; indicates whether the diagonal of the distance matrix should be included in the 
#' computation. Typically set to FALSE.
#' @param debugging Logical; set to FALSE by default. If TRUE, additional console messages will be 
#' displayed for debugging purposes.
#' @param isfloat Logical; indicates whether the input data values are floating-point numbers.
#' @param mfactor Integer; indicates the decimal position to round.
#' @return A list of matrices corresponding to the computed Rao's index values. Each matrix in the list 
#' represents the calculations performed over the moving window, with dimensions equal to \code{dim(x)}.
#' @author Duccio Rocchini \email{duccio.rocchini@@unibo.it},
#' Matteo Marcantonio \email{marcantoniomatteo@@gmail.com}
#' @seealso \code{\link{paRao}} for the related non-sequential function.

paRaoS <- function(x, alpha, window, dist_m, na.tolerance, diag, debugging, isfloat, mfactor) 
{
	# Some initial housekeeping
	  # `win` is the operative moving window
  win = window 
  NAwin <- 2*window+1
  message("\n\nProcessing alpha: ",alpha, " Moving Window: ", NAwin)
  # Set a progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent in :elapsed\n",
    # Total number of ticks is the number of column +NA columns divided the number of processor.
    total = dim(x)[2], 
    clear = FALSE, 
    width = 60, 
    force = FALSE)

	message("\n\nProcessing alpha: ",alpha, " Moving Window: ", 2*win+1)
	mfactor <- ifelse(isfloat,mfactor,1) 
	window = 2*win+1
	diagonal <- ifelse(diag==TRUE,0,NA)
	tdist <- proxy::dist(as.numeric(levels(as.factor(x))),method=dist_m)
	# Min and max dist for initial checks on possible infinite or 0 operations
	maxd <- max(proxy::dist(as.numeric(levels(as.factor(x))),method=dist_m))
	mind <- min(tdist[tdist>0])
	# If alpha ~ +infinite
	if( alpha >= .Machine$integer.max | is.infinite(alpha) | is.infinite(maxd^alpha) | (dist_m=="canberra" & mind^alpha==0) ) {
		paRaoOS <- matrix(rep(NA,dim(x)[1]*dim(x)[2]),nrow=dim(x)[1],ncol=dim(x)[2])
		# Reshape values
		values <- as.numeric(as.factor(x))
		x_1 <- matrix(data=values,nrow=dim(x)[1],ncol=dim(x)[2])
		# Add additional columns and rows for moving window
		hor <- matrix(NA,ncol=dim(x)[2],nrow=win)
		ver <- matrix(NA,ncol=win,nrow=dim(x)[1]+win*2)
		tx <- cbind(ver,rbind(hor,x_1,hor),ver)
		# Derive distance matrix
		classes <- levels(as.factor(x))
		if( is.character(dist_m) | is.function(dist_m) ) {
			d1 <- proxy::dist(as.numeric(classes),method=dist_m)
		}else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
			d1 <- stats::as.dist(stats::xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
		}
		# Loop over each pixel
		for (cl in (1+win):(dim(x)[2]+win)) {
			# Update progress bar
			pb$tick()
        	# Row loop
			for(rw in (1+win):(dim(x)[1]+win)) {
				if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) <= (window^2-((window^2)*na.tolerance)) ) {
					paRaoOS[rw-win,cl-win] <- NA
				}else{
					tw <- summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]),maxsum=10000)
					if( "NA's"%in%names(tw) ) {
						tw <- tw[-length(tw)]
					}
					if( debugging ) {
						message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window^2)
					}
					tw_labels <- names(tw)
					tw_values <- as.vector(tw)
          			# Exclude windows with only 1 category
					if( length(tw_values) == 1 ) {
						paRaoOS[rw-win,cl-win] <- 0
					}else{
						d2 <- unname(proxy::as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
						paRaoOS[rw-win,cl-win] <- max(d2*2,na.rm=TRUE) / mfactor
					}
				}
			}
		}
		return(paRaoOS)
	#If alpha is >0
	}else if( alpha > 0 ){
		paRaoOS <- matrix(rep(NA,dim(x)[1]*dim(x)[2]),nrow=dim(x)[1],ncol=dim(x)[2])
		# Reshape values
		values<-as.numeric(as.factor(x))
		x_1<-matrix(data=values,nrow=dim(x)[1],ncol=dim(x)[2])
		# Add additional columns and rows for moving window
		hor<-matrix(NA,ncol=dim(x)[2],nrow=win)
		ver<-matrix(NA,ncol=win,nrow=dim(x)[1]+win*2)
		tx<-cbind(ver,rbind(hor,x_1,hor),ver)
		# Derive distance matrix
		classes<-levels(as.factor(x))
		if( is.character(dist_m) | is.function(dist_m) ) {
			d1<-proxy::dist(as.numeric(classes),method=dist_m)
		}else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
			d1<-stats::as.dist(stats::xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
		}
		# Loop over each pixel
		for(cl in (1+win):(ncol(x)+win)) {
			# Update progress bar
			pb$tick()
			# Row loop
			for(rw in (1+win):(nrow(x)+win)) {
				if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) <= window^2-((window^2)*na.tolerance) ) {
					paRaoOS[rw-win,cl-win]<-NA
				}else{
					tw <- summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]),maxsum=10000)
					if( "NA's"%in%names(tw) ) {
						tw <- tw[-length(tw)]
					}
					if(debugging) {
						message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window^2)
					}
					tw_labels <- names(tw)
					tw_values <- as.vector(tw)
        			# if clause to exclude windows with less than 1 category
					if( length(tw_values) == 1 ) {
						paRaoOS[rw-win,cl-win] <- 0
					} else {
						p <- tw_values/sum(tw_values,na.rm=TRUE)
						p1 <- diag(diagonal,length(tw_values))
						p1[lower.tri(p1)] <- c(utils::combn(p,m=2,FUN=prod,na.rm=TRUE))
						d2 <- unname(proxy::as.matrix(d1)[as.numeric(tw_labels),as.numeric(tw_labels)])
						paRaoOS[rw-win,cl-win] <- (sum((p1)*(d2^alpha)*2,na.rm=TRUE))^(1/alpha) / mfactor
					}
				}
			}			
		} 
		return(paRaoOS)
	# If alpha == 0
	}else if( alpha==0 ){
		paRaoOS <- matrix(rep(NA,dim(x)[1]*dim(x)[2]),nrow=dim(x)[1],ncol=dim(x)[2])
			# Reshape values
		values <- as.numeric(as.factor(x))
		x_1 <- matrix(data=values,nrow=dim(x)[1],ncol=dim(x)[2])
			# Add additional columns and rows for moving window
		hor <- matrix(NA,ncol=dim(x)[2],nrow=win)
		ver <- matrix(NA,ncol=win,nrow=dim(x)[1]+win*2)
		tx <- cbind(ver,rbind(hor,x_1,hor),ver)
			# Derive distance matrix
		classes <- levels(as.factor(x))
		if( is.character(dist_m) | is.function(dist_m) ) {
			d1 <- proxy::dist(as.numeric(classes),method=dist_m)
		} else if( is.matrix(dist_m) | is.data.frame(dist_m) ) {
			d1 <- stats::as.dist(stats::xtabs(dist_m[, 3] ~ dist_m[, 2] + dist_m[, 1]))
		}
		for (cl in (1+win):(dim(x)[2]+win)) {
			for(rw in (1+win):(dim(x)[1]+win)) {
				if( length(!which(!tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]%in%NA)) <= (window^2-((window^2)*na.tolerance)) ) {
					paRaoOS[rw-win,cl-win] <- NA
				}else{
					tw <- summary(as.factor(tx[c(rw-win):c(rw+win),c(cl-win):c(cl+win)]),maxsum=10000)
					if( "NA's"%in%names(tw) ) {
						tw <- tw[-length(tw)]
					}
					if(debugging) {
						message("Working on coords ",rw ,",",cl,". classes length: ",length(tw),". window size=",window^2)
					}
					tw_labels <- names(tw)
					tw_values <- as.vector(tw)
      					# if clause to exclude windows with only 1 category
					if(length(tw_values) == 1) {
						paRaoOS[rw-win,cl-win] <- 0
					}else{
						d2 <- unname( proxy::as.matrix(d1,diag=diagonal)[as.numeric(tw_labels),as.numeric(tw_labels)] )
						paRaoOS[rw-win,cl-win] <- ( prod(d2/mfactor,na.rm=TRUE)^(1/(window^2)) )
					}
				} 
			} 
		}
		return(paRaoOS)
	}
}