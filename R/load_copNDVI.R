#' Load Copernicus NDVI data
#'
#' This function loads and returns the Copernicus Long Term (1999-2017) NDVI Overview stored within the package.
#'
#' @return A `SpatRaster` object representing the Copernicus NDVI data.
#' @export
#' @examples
#' copNDVI <- load_copNDVI()
load_copNDVI <- function() {
  # Set the path to the RDS file
  rds_file <- system.file("extdata", "copNDVI.rds", package = "rasterdiv")
  
  if (file.exists(rds_file)) {
    copNDVI <- terra::rast(readRDS(rds_file))
  } else {
    stop("Data file 'copNDVI.rds' not found in the package.")
  }
  
  return(copNDVI)
}
