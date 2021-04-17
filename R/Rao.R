Rao = function(x, dist_m="euclidean", window=9, rasterOut = TRUE, mode="classic",lambda=0, shannon=FALSE, rescale=FALSE, na.tolerance=1.0, simplify=2, np=1, cluster.type="SOCK",debugging=FALSE) {
	.Deprecated(new = "paRao(..., alpha=1)")
	paRao(x=x, dist_m=dist_m, window=window, method=mode, alpha=1, lambda=lambda, na.tolerance=na.tolerance, rescale=rescale, diag=TRUE, simplify=simplify, np=np, cluster.type=cluster.type, debugging=debugging)
}