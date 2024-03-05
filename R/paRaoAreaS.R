#' Area-Based Sequential Parametric Rao's index of quadratic entropy (Q)
#'
#' Area-Based Sequential Parametric Rao's index of quadratic entropy (Q).
#'
#' @param rasterm input raster data.
#' @param area input vector data.
#' @param alpha alpha value for order of diversity in Hill's Index.
#' @param simplify Rounding parameter.
#'
#' @return A vector as the input with Rao's index pasted additional columns.
#'
#' @author Marcantonio Matteo \email{marcantoniomatteo@gmail.com}, 
#' Duccio Rocchini \email{duccio.rocchini@unibo.it}, 
#' Michele Torresani \email{michele.torresani@unibo.it}
#'
#' @seealso \code{\link{paRao}}
#' 
#' @keywords internal

paRaoAreaS <- function(rasterm, area, alpha, simplify) {
	if(alpha<=0) {
		stop("Alpha<=0 not yet implemented for areas.")
	}
	if(is.matrix(rasterm)) {
		stop("x must be a SpatRaster.")
	}

	mfactor <- 100^simplify
	crop1 <- terra::crop(rasterm, area)
	crop1dt <- terra::as.matrix(crop1)*mfactor
	storage.mode(crop1dt) <- "integer"

# Check for only 1 value in the matrix
if( any(is.na(crop1dt)) & all(apply(crop1dt, 2, function(a) length(unique(a))<=2)) ) {
	paRaoOareaS <- NA
	} else {
		classes <- levels(as.factor(crop1dt))
		d1 <- unname(proxy::as.matrix(proxy::dist(as.numeric(classes),method="Euclidean")))
		tw <- summary(as.factor(crop1dt),maxsum=10000)
		if( "NA's"%in%names(tw) ) {
			tw <- tw[-length(tw)]
		}

		tw_values <- as.vector(tw)
		p <- tw_values/sum(tw_values, na.rm=TRUE)
		p1 <- diag(TRUE,length(tw_values))
		p1[lower.tri(p1)] <- c(utils::combn(p,m=2,FUN=prod,na.rm=TRUE))
		paRaoOareaS <- (sum((p1)*(d1^alpha)*2,na.rm=TRUE))^(1/alpha) / mfactor
		gc()
	}
	return(paRaoOareaS)
}