paRaoAreaS <- function(rasterm, area, alpha) {
 if(alpha<=0) {stop("Alpha<=0 not yet implemented for areas")}
 if(is.matrix(rasterm)) {stop("Rasterm must be a raster")}
  crop1 <- crop(rasterm, area)
  mat_s <- values(crop1)
  n_s <- length(mat_s)
  n2_s <- n_s^alpha
  distm_s <- as.matrix(dist(mat_s))
  rao <- (sum(distm_s, na.rm=TRUE))/n2_s
  gc()
  return(rao)
}