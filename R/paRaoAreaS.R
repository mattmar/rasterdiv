paRaoAreaS <- function(rasterm, area, alpha, simplify) {
	if(alpha<=0) {
		stop("Alpha<=0 not yet implemented for areas.")
	}
	if(is.matrix(rasterm)) {
		stop("x must be a RasterLayer.")
	}

	mfactor <- 100^simplify
	crop1 <- crop(rasterm, area)
	crop1dt <- raster::as.matrix(crop1)*mfactor
	storage.mode(crop1dt) <- "integer"
	classes <- levels(as.factor(crop1dt))
	d1 <- unname(proxy::as.matrix(proxy::dist(as.numeric(classes),method="Euclidean")))

	tw <- summary(as.factor(crop1dt),maxsum=10000)
	if( "NA's"%in%names(tw) ) {
		tw <- tw[-length(tw)]
	}

	tw_values <- as.vector(tw)
	p <- tw_values/sum(tw_values, na.rm=TRUE)
	p1 <- diag(TRUE,length(tw_values))
	p1[lower.tri(p1)] <- c(combn(p,m=2,FUN=prod,na.rm=TRUE))
	paRaoOareaS <- (sum((p1)*(d1^alpha)*2,na.rm=TRUE))^(1/alpha) / mfactor
	gc()
	return(paRaoOareaS)
}