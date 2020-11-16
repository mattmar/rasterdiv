accRao <- function(alphas=1:5, x, dist_m="euclidean", window=9, method="classic", rasterAUC=TRUE, lambda=0, na.tolerance=1.0, rescale=FALSE, diag=TRUE, simplify=2, np=1, cluster.type="SOCK", debugging=FALSE)
{
	out <- paRao(x, dist_m, window, method, rasterOut=FALSE, alpha=alphas, lambda, na.tolerance, rescale, diag, simplify, np, cluster.type, debugging)

	message("\nIntegrating numerically Rao values over alphas...\n")

	outafx <- lapply(out, function(w) {
		outm <- apply(
			sapply(w, function(a) {
				sapply(1:length(a), function(y) {
					a[y]
				})
			})
			,1, function(i) {
				if( all(is.na(i)) ) {
					NA
				} else {
					integrate(approxfun(y=i,x=alphas), lower=alphas[1], upper=alphas[max(alphas)], subdivisions = 500)$value
				}
			})
		return(outm)
	})

	if( rasterAUC==TRUE & class(x)[[1]]=="RasterLayer" ) {
		outR <- lapply(outafx, function(insm) {
			y <- raster(matrix(insm,ncol=ncol(x),nrow=nrow(x)),template=x)
		})
		return(outR)
	} else {
		outM <- lapply(outafx, function(insm) {
			y <- matrix(insm,ncol=ncol(x),nrow=nrow(x))
		})
		return(outM)
	}
}